/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FlowSolverBase.cpp
 */

#include "FlowSolverBase.hpp"

#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/solid/SolidInternalEnergy.hpp"
#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/LogLevelsInfo.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/MinPoreVolumeMaxPorosityKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/StencilWeightsUpdateKernel.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

template< typename POROUSWRAPPER_TYPE >
void updatePorosityAndPermeabilityFromPressureAndTemperature( POROUSWRAPPER_TYPE porousWrapper,
                                                              CellElementSubRegion & subRegion,
                                                              arrayView1d< real64 const > const & pressure,
                                                              arrayView1d< real64 const > const & temperature )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFromPressureAndTemperature( k, q,
                                                           pressure[k],
                                                           temperature[k] );
    }
  } );
}

template< typename POROUSWRAPPER_TYPE >
void updatePorosityAndPermeabilityFixedStress( POROUSWRAPPER_TYPE porousWrapper,
                                               CellElementSubRegion & subRegion,
                                               arrayView1d< real64 const > const & pressure,
                                               arrayView1d< real64 const > const & pressure_k,
                                               arrayView1d< real64 const > const & pressure_n,
                                               arrayView1d< real64 const > const & temperature,
                                               arrayView1d< real64 const > const & temperature_k,
                                               arrayView1d< real64 const > const & temperature_n )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFixedStress( k, q,
                                            pressure[k],
                                            pressure_k[k],
                                            pressure_n[k],
                                            temperature[k],
                                            temperature_k[k],
                                            temperature_n[k] );
    }
  } );
}

template< typename POROUSWRAPPER_TYPE >
void updatePorosityAndPermeabilityFromPressureAndAperture( POROUSWRAPPER_TYPE porousWrapper,
                                                           SurfaceElementSubRegion & subRegion,
                                                           arrayView1d< real64 const > const & pressure,
                                                           arrayView1d< real64 const > const & oldHydraulicAperture,
                                                           arrayView1d< real64 const > const & newHydraulicAperture )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFromPressureAndAperture( k, q,
                                                        pressure[k],
                                                        oldHydraulicAperture[k],
                                                        newHydraulicAperture[k] );
    }
  } );
}

FlowSolverBase::FlowSolverBase( string const & name,
                                Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_numDofPerCell( 0 ),
  m_isThermal( 0 ),
  m_keepVariablesConstantDuringInitStep( false ),
  m_isFixedStressPoromechanicsUpdate( false ),
  m_isJumpStabilized( false ),
  m_isLaggingFractureStencilWeightsUpdate( 0 )
{
  this->registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not." );

  this->registerWrapper( viewKeyStruct::allowNegativePressureString(), &m_allowNegativePressure ).
    setApplyDefaultValue( 0 ). // negative pressure is not allowed by default
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating if negative pressure is allowed" );

  this->registerWrapper( viewKeyStruct::maxAbsolutePresChangeString(), &m_maxAbsolutePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).       // disabled by default
    setDescription( "Maximum (absolute) pressure change in a Newton iteration" );

  this->registerWrapper( viewKeyStruct::maxSequentialPresChangeString(), &m_maxSequentialPresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1e5 ).           // 0.1 bar = 1e5 Pa
    setDescription( "Maximum (absolute) pressure change in a sequential iteration, used for outer loop convergence check" );

  this->registerWrapper( viewKeyStruct::maxSequentialTempChangeString(), &m_maxSequentialTempChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription( "Maximum (absolute) temperature change in a sequential iteration, used for outer loop convergence check" );

  // allow the user to select a norm
  getNonlinearSolverParameters().getWrapper< physicsSolverBaseKernels::NormType >( NonlinearSolverParameters::viewKeysStruct::normTypeString() ).setInputFlag( InputFlags::OPTIONAL );

  addLogLevel< logInfo::Convergence >();
}

void FlowSolverBase::registerDataOnMesh( Group & meshBodies )
{
  PhysicsSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< fields::flow::deltaVolume >( getName() );
      subRegion.registerField< fields::flow::gravityCoefficient >( getName() ).
        setApplyDefaultValue( 0.0 );
      subRegion.registerField< fields::flow::netToGross >( getName() );

      subRegion.registerField< fields::flow::pressure >( getName() );
      subRegion.registerField< fields::flow::pressure_n >( getName() );
      subRegion.registerField< fields::flow::initialPressure >( getName() );
      subRegion.registerField< fields::flow::deltaPressure >( getName() ); // for reporting/stats purposes
      subRegion.registerField< fields::flow::bcPressure >( getName() ); // needed for the application of boundary conditions
      if( m_isFixedStressPoromechanicsUpdate )
      {
        subRegion.registerField< fields::flow::pressure_k >( getName() ); // needed for the fixed-stress porosity update
      }

      subRegion.registerField< fields::flow::temperature >( getName() );
      subRegion.registerField< fields::flow::temperature_n >( getName() );
      subRegion.registerField< fields::flow::initialTemperature >( getName() );
      subRegion.registerField< fields::flow::bcTemperature >( getName() ); // needed for the application of boundary conditions
      if( m_isFixedStressPoromechanicsUpdate )
      {
        subRegion.registerField< fields::flow::temperature_k >( getName() ); // needed for the fixed-stress porosity update
      }
      if( m_isThermal )
      {
        subRegion.registerField< fields::flow::energy >( getName() );
        subRegion.registerField< fields::flow::energy_n >( getName() );
      }
    } );

    elemManager.forElementSubRegionsComplete< SurfaceElementSubRegion >( [&]( localIndex const,
                                                                              localIndex const,
                                                                              ElementRegionBase & region,
                                                                              SurfaceElementSubRegion & subRegion )
    {
      SurfaceElementRegion & faceRegion = dynamicCast< SurfaceElementRegion & >( region );

      subRegion.registerField< fields::flow::gravityCoefficient >( getName() );

      subRegion.registerField< fields::flow::aperture0 >( getName() ).
        setApplyDefaultValue( faceRegion.getDefaultAperture() );

      subRegion.registerField< fields::flow::hydraulicAperture >( getName() ).
        setApplyDefaultValue( faceRegion.getDefaultAperture() );

    } );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< fields::flow::gravityCoefficient >( getName() ).
      setApplyDefaultValue( 0.0 );
    faceManager.registerField< fields::flow::transMultiplier >( getName() );

  } );

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {

    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
    fluxApprox.addFieldName( fields::flow::pressure::key() );
    fluxApprox.setCoeffName( fields::permeability::permeability::key() );
    if( m_isThermal )
    {
      fluxApprox.addFieldName( fields::flow::temperature::key() );
    }
  }
}

void FlowSolverBase::saveConvergedState( ElementSubRegionBase & subRegion ) const
{
  arrayView1d< real64 const > const pres = subRegion.template getField< fields::flow::pressure >();
  arrayView1d< real64 > const pres_n = subRegion.template getField< fields::flow::pressure_n >();
  pres_n.setValues< parallelDevicePolicy<> >( pres );

  arrayView1d< real64 const > const temp = subRegion.template getField< fields::flow::temperature >();
  arrayView1d< real64 > const temp_n = subRegion.template getField< fields::flow::temperature_n >();
  temp_n.setValues< parallelDevicePolicy<> >( temp );

  if( m_isThermal )
  {
    arrayView1d< real64 const > const energy = subRegion.template getField< fields::flow::energy >();
    arrayView1d< real64 > const energy_n = subRegion.template getField< fields::flow::energy_n >();
    energy_n.setValues< parallelDevicePolicy<> >( energy );
  }

  if( m_isFixedStressPoromechanicsUpdate )
  {
    arrayView1d< real64 > const pres_k = subRegion.template getField< fields::flow::pressure_k >();
    arrayView1d< real64 > const temp_k = subRegion.template getField< fields::flow::temperature_k >();
    pres_k.setValues< parallelDevicePolicy<> >( pres );
    temp_k.setValues< parallelDevicePolicy<> >( temp );
  }
}

void FlowSolverBase::saveSequentialIterationState( DomainPartition & domain )
{
  GEOS_ASSERT( m_isFixedStressPoromechanicsUpdate );

  real64 maxPresChange = 0.0;
  real64 maxTempChange = 0.0;
  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&]( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions ( regionNames,
                                                 [&]( localIndex const,
                                                      ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 > const pres_k = subRegion.getField< fields::flow::pressure_k >();
      arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
      arrayView1d< real64 > const temp_k = subRegion.getField< fields::flow::temperature_k >();

      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxPresChange( 0.0 );
      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxTempChange( 0.0 );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          subRegionMaxPresChange.max( LvArray::math::abs( pres[ei] - pres_k[ei] ) );
          pres_k[ei] = pres[ei];
          subRegionMaxTempChange.max( LvArray::math::abs( temp[ei] - temp_k[ei] ) );
          temp_k[ei] = temp[ei];
        }
      } );

      maxPresChange = LvArray::math::max( maxPresChange, subRegionMaxPresChange.get() );
      maxTempChange = LvArray::math::max( maxTempChange, subRegionMaxTempChange.get() );
    } );
  } );

  // store to be later used in convergence check
  m_sequentialPresChange = MpiWrapper::max( maxPresChange );
  m_sequentialTempChange = m_isThermal ? MpiWrapper::max( maxTempChange ) : 0.0;
}

void FlowSolverBase::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  PhysicsSolverBase::setConstitutiveNamesCallSuper( subRegion );

  setConstitutiveName< CoupledSolidBase >( subRegion, viewKeyStruct::solidNamesString(), "coupled solid" );

  setConstitutiveName< PermeabilityBase >( subRegion, viewKeyStruct::permeabilityNamesString(), "permeability" );

  if( m_isThermal )
  {
    setConstitutiveName< SolidInternalEnergy >( subRegion, viewKeyStruct::solidInternalEnergyNamesString(), "solid internal energy" );
  }
}

void FlowSolverBase::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  GEOS_UNUSED_VAR( subRegion );
}

void FlowSolverBase::initializePreSubGroups()
{
  PhysicsSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                  MeshLevel &,
                                                                  string_array const & regionNames )
    {
      string_array & stencilTargetRegions = fluxApprox.targetRegions( meshBodyName );
      std::set< string > stencilTargetRegionsSet( stencilTargetRegions.begin(), stencilTargetRegions.end() );
      stencilTargetRegionsSet.insert( regionNames.begin(), regionNames.end() );
      stencilTargetRegions.clear();
      for( auto const & targetRegion: stencilTargetRegionsSet )
      {
        stencilTargetRegions.emplace_back( targetRegion );
      }
    } );
  }
}

void FlowSolverBase::checkDiscretizationName() const
{
  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOS_ERROR( GEOS_FMT( "{}: can not find discretization named '{}' (a discretization deriving from FluxApproximationBase must be selected for {} solver '{}' )",
                          getDataContext(), m_discretizationName, getCatalogName(), getName()));
  }
}

void FlowSolverBase::validatePoreVolumes( DomainPartition const & domain ) const
{
  real64 minPoreVolume = LvArray::NumericLimits< real64 >::max;
  real64 maxPorosity = -LvArray::NumericLimits< real64 >::max;
  globalIndex numElemsBelowPoreVolumeThreshold = 0;
  globalIndex numElemsAbovePorosityThreshold = 0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase const & subRegion )
    {

      string const & solidName = subRegion.template getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      real64 minPoreVolumeInSubRegion = 0.0;
      real64 maxPorosityInSubRegion = 0.0;
      localIndex numElemsBelowPoreVolumeThresholdInSubRegion = 0;
      localIndex numElemsAbovePorosityThresholdInSubRegion = 0;

      flowSolverBaseKernels::MinPoreVolumeMaxPorosityKernel::
        computeMinPoreVolumeMaxPorosity( subRegion.size(),
                                         ghostRank,
                                         porosity,
                                         volume,
                                         minPoreVolumeInSubRegion,
                                         maxPorosityInSubRegion,
                                         numElemsBelowPoreVolumeThresholdInSubRegion,
                                         numElemsAbovePorosityThresholdInSubRegion );

      if( minPoreVolumeInSubRegion < minPoreVolume )
      {
        minPoreVolume = minPoreVolumeInSubRegion;
      }
      if( maxPorosityInSubRegion > maxPorosity )
      {
        maxPorosity = maxPorosityInSubRegion;
      }

      numElemsBelowPoreVolumeThreshold += numElemsBelowPoreVolumeThresholdInSubRegion;
      numElemsAbovePorosityThreshold += numElemsAbovePorosityThresholdInSubRegion;
    } );
  } );

  minPoreVolume = MpiWrapper::min( minPoreVolume );
  maxPorosity = MpiWrapper::max( maxPorosity );
  numElemsBelowPoreVolumeThreshold = MpiWrapper::sum( numElemsBelowPoreVolumeThreshold );
  numElemsAbovePorosityThreshold = MpiWrapper::sum( numElemsAbovePorosityThreshold );

  GEOS_LOG_RANK_0_IF( numElemsBelowPoreVolumeThreshold > 0,
                      GEOS_FMT( "\nWarning! The mesh contains {} elements with a pore volume below {} m^3."
                                "\nThe minimum pore volume is {} m^3."
                                "\nOur recommendation is to check the validity of mesh and/or increase the porosity in these elements.\n",
                                numElemsBelowPoreVolumeThreshold, flowSolverBaseKernels::poreVolumeThreshold, minPoreVolume ) );
  GEOS_LOG_RANK_0_IF( numElemsAbovePorosityThreshold > 0,
                      GEOS_FMT( "\nWarning! The mesh contains {} elements with a porosity above 1."
                                "\nThe maximum porosity is {}.\n",
                                numElemsAbovePorosityThreshold, maxPorosity ) );
}

void FlowSolverBase::initializePostInitialConditionsPreSubGroups()
{
  PhysicsSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    precomputeData( mesh, regionNames );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::flow::pressure::key(), fields::flow::temperature::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );
  } );
}

void FlowSolverBase::precomputeData( MeshLevel & mesh,
                                     string_array const & regionNames )
{
  FaceManager & faceManager = mesh.getFaceManager();
  real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                        ElementSubRegionBase & subRegion )
  {
    arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();

    arrayView1d< real64 > const gravityCoef =
      subRegion.getField< fields::flow::gravityCoefficient >();

    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      gravityCoef[ ei ] = LvArray::tensorOps::AiBi< 3 >( elemCenter[ ei ], gravVector );
    } );
  } );

  {
    arrayView2d< real64 const > const faceCenter = faceManager.faceCenter();

    arrayView1d< real64 > const gravityCoef =
      faceManager.getField< fields::flow::gravityCoefficient >();

    forAll< parallelHostPolicy >( faceManager.size(), [=] ( localIndex const kf )
    {
      gravityCoef[ kf ] = LvArray::tensorOps::AiBi< 3 >( faceCenter[ kf ], gravVector );
    } );
  }
}

void FlowSolverBase::initializeState( DomainPartition & domain )
{
  // Compute hydrostatic equilibrium in the regions for which corresponding field specification tag has been specified
  computeHydrostaticEquilibrium( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    initializePorosityAndPermeability( mesh, regionNames );
    initializeHydraulicAperture( mesh, regionNames );

    // Initialize primary variables from applied initial conditions
    initializeFluidState( mesh, regionNames );

    // Initialize the rock thermal quantities: conductivity and solid internal energy
    // Note:
    // - This must be called after updatePorosityAndPermeability and updatePhaseVolumeFraction
    // - This step depends on porosity and phaseVolFraction
    if( m_isThermal )
    {
      initializeThermalState( mesh, regionNames );
    }

    // Save initial pressure and temperature fields
    saveInitialPressureAndTemperature( mesh, regionNames );
  } );

  // report to the user if some pore volumes are very small
  // note: this function is here because: 1) porosity has been initialized and 2) NTG has been applied
  validatePoreVolumes( domain );
}

void FlowSolverBase::initializePorosityAndPermeability( MeshLevel & mesh, string_array const & regionNames )
{
  // Update porosity and permeability
  // In addition, to avoid multiplying permeability/porosity bay netToGross in the assembly kernel, we do it once and for all here
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                 auto & subRegion )
  {
    // Apply netToGross to reference porosity and horizontal permeability
    CoupledSolidBase const & porousSolid =
      getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
    PermeabilityBase const & permeability =
      getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() ) );
    arrayView1d< real64 const > const netToGross = subRegion.template getField< fields::flow::netToGross >();
    porousSolid.scaleReferencePorosity( netToGross );
    permeability.scaleHorizontalPermeability( netToGross );

    // in some initializeState versions it uses newPorosity, so let's run updatePorosityAndPermeability to compute something
    saveConvergedState( subRegion );   // necessary for a meaningful porosity update in sequential schemes
    updatePorosityAndPermeability( subRegion );
    porousSolid.initializeState();

    // run final update
    saveConvergedState( subRegion );   // necessary for a meaningful porosity update in sequential schemes
    updatePorosityAndPermeability( subRegion );

    // Save the computed porosity into the old porosity
    // Note:
    // - This must be called after updatePorosityAndPermeability
    // - This step depends on porosity
    porousSolid.saveConvergedState();
  } );
}

void FlowSolverBase::initializeHydraulicAperture( MeshLevel & mesh, string_array const & regionNames )
{
  mesh.getElemManager().forElementRegions< SurfaceElementRegion >( regionNames,
                                                                   [&]( localIndex const,
                                                                        SurfaceElementRegion & region )
  {
    region.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    { subRegion.getWrapper< real64_array >( fields::flow::hydraulicAperture::key()).setApplyDefaultValue( region.getDefaultAperture()); } );
  } );
}

void FlowSolverBase::saveInitialPressureAndTemperature( MeshLevel & mesh, string_array const & regionNames )
{
  mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 > const initPres = subRegion.getField< fields::flow::initialPressure >();
    arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
    arrayView1d< real64 > const initTemp = subRegion.template getField< fields::flow::initialTemperature >();
    initPres.setValues< parallelDevicePolicy<> >( pres );
    initTemp.setValues< parallelDevicePolicy<> >( temp );
  } );
}

void FlowSolverBase::updatePorosityAndPermeability( CellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const & pressure = subRegion.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const & temperature = subRegion.getField< fields::flow::temperature >();
  // LILIANE  
  /*forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  {
    GEOS_LOG_RANK_0( "*** "<< k <<" matrix pressure "<< pressure[k] );
  });*/
  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( solidName );

  constitutive::ConstitutivePassThru< CoupledSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
  {
    typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();
    if( m_isFixedStressPoromechanicsUpdate )
    {
      arrayView1d< real64 const > const & pressure_n = subRegion.getField< fields::flow::pressure_n >();
      arrayView1d< real64 const > const & pressure_k = subRegion.getField< fields::flow::pressure_k >();
      arrayView1d< real64 const > const & temperature_n = subRegion.getField< fields::flow::temperature_n >();
      arrayView1d< real64 const > const & temperature_k = subRegion.getField< fields::flow::temperature_k >();
      updatePorosityAndPermeabilityFixedStress( porousWrapper, subRegion, pressure, pressure_k, pressure_n, temperature, temperature_k, temperature_n );
    }
    else
    {
      updatePorosityAndPermeabilityFromPressureAndTemperature( porousWrapper, subRegion, pressure, temperature );
    }
  } );
}

// LILIANE
void updateFractureAperture_old( SurfaceElementSubRegion & subRegion ) 
{
  GEOS_MARK_FUNCTION;

  auto dot_product = [](real64 const (&u)[3], real64 const (&v)[3] ) -> real64
  {
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
  };
  auto matmul = [](real64 const (&u)[3], real64 const (&v)[3], real64 (&r)[3]) -> void
  {
    r[0] = u[0]*v[0];
    r[1] = u[1]*v[1];
    r[2] = u[2]*v[2];
  };

  arrayView1d< real64 const > const & pressure = subRegion.getField< fields::flow::pressure >();
  /*if (pressure.size() > 0)
  {
    forAll< parallelDevicePolicy<> >( subRegion.size(), [&] GEOS_DEVICE ( localIndex const k )
    {
      GEOS_LOG_RANK_0( "*** "<< k << std::setprecision(6) << std::scientific <<" pressure "<< pressure[k]);
    } );
  }*/
  arrayView1d< real64 > newHydraulicAperture = subRegion.getField< fields::flow::hydraulicAperture >();
  auto normalVector = subRegion.getField< fields::normalVector >();
  // Estamos considerando que o estado de referência é o estado in-situ (reservatório)

  // As variáveis abaixo serão lidas, num próximo passo, do xml:
  // Tensão total no estado de referência
  real64 const sigmaT_0[3] = {85.0e6, 85.0e6, 105.0e6 };     

  // pressão no estado de referência
  // esta deve ser a condição inicial de pressão do problema a ser resolvido (ajustar o xml do problema que 
  // você  está rodando agora
  //real64 const p_0  = 55.0e6; 
  real64 const p_0  = 1e5;             

  real64 const biot = 1.0;       // coeficiente de Biot
  real64 const poisson = 0.3;    // coeficiente de poisson
  real64 const Kni = 12.041e9;   // rigidez normal inicial

  // Vm é o fechamento máximo, igual ao  valor da abertura no estado livre de tensões (pode variar por fratura)
  real64 const Vm = 1e-3; 
  
  // Cálculos a serem feitos para calcular a abertura da fratura para uma dada pressão 
  // Tem que ser feitos para cada elemento de cada fratura (inicialmente vamos considerar uma única fratura)

  // sigma_0 é a tensão efetiva no estado de referência
  real64 const sigma_0[3] = {sigmaT_0[0]-biot*p_0, sigmaT_0[1]-biot*p_0, sigmaT_0[2]-biot*p_0};
  // sigma é a tensão efetiva na pressão dada; varia de acordo com a trajetória edométrica que assumimos aqui

  real64 sumAperture = 0.0;
  real64 sumSigmaN = 0.0;
  forAll< parallelDevicePolicy<> >( subRegion.size(), [&] GEOS_DEVICE ( localIndex const k )
  {
    //for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      real64 const delta_sigmaZ = biot*(pressure[k]-p_0);     // pressao_Atual é a pressão no elemento da fratura
      real64 sigma[3] = { sigma_0[0] - (poisson/(1.0-poisson)) * delta_sigmaZ,
                          sigma_0[1] - (poisson/(1.0-poisson)) * delta_sigmaZ,
                          sigma_0[2] - delta_sigmaZ };                          
      // Calcular a normal à fratura 
      // Calcular a tensão efetiva que vai atuar sobre a fratura
      
      real64 normalElem[3] = {normalVector[k][0], normalVector[k][1], normalVector[k][2]};
      matmul(sigma,normalElem, sigma);
      real64 const sigmaN_N = dot_product(sigma, normalElem);
      
      // Calcular o fechamento através da lei constitutiva de BB
      real64 const gn_BB =  sigmaN_N*Vm/(Kni*Vm + sigmaN_N);

      // Calcular a nova abertura que é igual à abertura no estado livre de tensões menos o fechemento a partir do 
      // estado livre de tensões até o estado atual
      newHydraulicAperture[k] = Vm - gn_BB;

      sumAperture += newHydraulicAperture[k];
	    sumSigmaN += sigmaN_N;
      
      GEOS_LOG_RANK_0( "*** "<< k << std::setprecision(6) << std::scientific <<" pressure "<< pressure[k] << " newHydraulicAperture "<< newHydraulicAperture[k]<< " sigmaN_N "<< sigmaN_N );
    }
  } );                      

  //No final do loop de cada fratura, calcular a média da abertura e atribuir à abertura de cada elemento daquela fratura
  real64 const averageAperture = sumAperture / subRegion.size();
  real64 const averageSigmaN = sumSigmaN / subRegion.size();
  forAll< parallelDevicePolicy<> >( subRegion.size(), [&newHydraulicAperture, averageAperture] GEOS_DEVICE ( localIndex const k )
  {
    newHydraulicAperture[k] = averageAperture;
  } ); 

  GEOS_LOG_RANK_0( std::setprecision(6) << std::scientific <<"*** averageAperture "<< averageAperture << " averageSigmaN "<< averageSigmaN );
  
}


void FlowSolverBase::updateFractureAperture( SurfaceElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  auto dot_product = [](real64 const (&u)[3], real64 const (&v)[3] ) -> real64
  {
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
  };
  auto matmul = [](real64 const (&u)[3], real64 const (&v)[3], real64 (&r)[3]) -> void
  {
    r[0] = u[0]*v[0];
    r[1] = u[1]*v[1];
    r[2] = u[2]*v[2];
  };

  arrayView1d< real64 const > const & pressure = subRegion.getField< fields::flow::pressure >();
  /*if (pressure.size() > 0)
  {
    forAll< parallelDevicePolicy<> >( subRegion.size(), [&] GEOS_DEVICE ( localIndex const k )
    {
      GEOS_LOG_RANK_0( "*** "<< k << std::setprecision(6) << std::scientific <<" pressure "<< pressure[k]);
    } );
  }*/
  arrayView1d< real64 > newHydraulicAperture = subRegion.getField< fields::flow::hydraulicAperture >();
  auto normalVector = subRegion.getField< fields::normalVector >();
  // Estamos considerando que o estado de referência é o estado in-situ (reservatório)

  // As variáveis abaixo serão lidas, num próximo passo, do xml:
  // Tensão total no estado de referência
  real64 const sigmaT_0[3] = {85.0e6, 85.0e6, 105.0e6 };     

  // pressão no estado de referência
  // esta deve ser a condição inicial de pressão do problema a ser resolvido (ajustar o xml do problema que 
  // você  está rodando agora
  //real64 const p_0  = 55.0e6; 
  real64 const p_0  = 1e5;             

  real64 const biot = 1.0;       // coeficiente de Biot
  real64 const poisson = 0.3;    // coeficiente de poisson
  real64 const Kni = 12.041e9;   // rigidez normal inicial

  // Vm é o fechamento máximo, igual ao  valor da abertura no estado livre de tensões (pode variar por fratura)
  //real64 const Vm = 1e-3; 

  // d0 é a abertura da fratura no estado in-situ
  real64 const d0 = 1.e-3; 
  GEOS_LOG_RANK_0( "d0 "<< k << std::setprecision(6) << std::scientific << d0);
  
  
  // Cálculos a serem feitos para calcular a abertura da fratura para uma dada pressão 
  // Tem que ser feitos para cada elemento de cada fratura (inicialmente vamos considerar uma única fratura)

  // sigma_0 é a tensão efetiva no estado de referência
  real64 const sigma_0[3] = {sigmaT_0[0]-biot*p_0, sigmaT_0[1]-biot*p_0, sigmaT_0[2]-biot*p_0};
  // sigma é a tensão efetiva na pressão dada; varia de acordo com a trajetória edométrica que assumimos aqui

  real64 sumAperture = 0.0;
  real64 sumSigmaN = 0.0;
  forAll< parallelDevicePolicy<> >( subRegion.size(), [&] GEOS_DEVICE ( localIndex const k )
  {
    //for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      real64 const delta_sigmaZ = biot*(pressure[k]-p_0);     // pressao_Atual é a pressão no elemento da fratura
      real64 sigma[3] = { sigma_0[0] - (poisson/(1.0-poisson)) * delta_sigmaZ,
                          sigma_0[1] - (poisson/(1.0-poisson)) * delta_sigmaZ,
                          sigma_0[2] - delta_sigmaZ };

      real64 const normalElem[3] = {normalVector[k][0], normalVector[k][1], normalVector[k][2]};
      real64 const alpha_p = biot * p_0;
      real64 const sigma_c[3] = { sigmaT_0[0] * normalElem[0] - alpha_p * normalElem[0],
                                  sigmaT_0[1] * normalElem[1] - alpha_p * normalElem[1],
                                  sigmaT_0[2] * normalElem[2] - alpha_p * normalElem[2]};
                                  
      real64 const sigma_n0 = sigma_c[0]*normalElem[0] + sigma_c[1]*normalElem[1] + sigma_c[2]*normalElem[2];
      real64 const g0 = (-Kni*d0 + std::sqrt((Kni*d0)*(Kni*d0)+4.0*Kni*sigma_n0*d0))/(2.0*Kni);
      real64 const Vm = g0 + d0;

                          
      // Calcular a normal à fratura 
      // Calcular a tensão efetiva que vai atuar sobre a fratura
      matmul(sigma, normalElem, sigma);
      real64 const sigmaN_N = dot_product(sigma, normalElem);
      
      // Calcular o fechamento através da lei constitutiva de BB
      real64 const gn_BB =  sigmaN_N*Vm/(Kni*Vm + sigmaN_N);

      // Calcular a nova abertura que é igual à abertura no estado livre de tensões menos o fechemento a partir do 
      // estado livre de tensões até o estado atual
      newHydraulicAperture[k] = Vm - gn_BB;

      sumAperture += newHydraulicAperture[k];
	    sumSigmaN += sigmaN_N;
      
      GEOS_LOG_RANK_0( "*** "<< k << std::setprecision(6) << std::scientific <<" pressure "<< pressure[k] << " newHydraulicAperture "<< newHydraulicAperture[k]<< " sigmaN_N "<< sigmaN_N << " sigma_n0 " << sigma_n0 << " Vm " << Vm);
    }
  } );                      

  //No final do loop de cada fratura, calcular a média da abertura e atribuir à abertura de cada elemento daquela fratura
  real64 const averageAperture = sumAperture / subRegion.size();
  real64 const averageSigmaN = sumSigmaN / subRegion.size();
  forAll< parallelDevicePolicy<> >( subRegion.size(), [&newHydraulicAperture, averageAperture] GEOS_DEVICE ( localIndex const k )
  {
    newHydraulicAperture[k] = averageAperture;
  } ); 

  GEOS_LOG_RANK_0( std::setprecision(6) << std::scientific <<"*** averageAperture "<< averageAperture << " averageSigmaN "<< averageSigmaN );
  
}

void FlowSolverBase::updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const & pressure = subRegion.getField< fields::flow::pressure >();

  arrayView1d< real64 const > const newHydraulicAperture = subRegion.getField< fields::flow::hydraulicAperture >();
  arrayView1d< real64 const > const oldHydraulicAperture = subRegion.getField< fields::flow::aperture0 >();

  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  CoupledSolidBase & porousSolid = subRegion.getConstitutiveModel< CoupledSolidBase >( solidName );

  constitutive::ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
  {
    typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();

    updatePorosityAndPermeabilityFromPressureAndAperture( porousWrapper, subRegion, pressure, oldHydraulicAperture, newHydraulicAperture );

  } );
}


void FlowSolverBase::findMinMaxElevationInEquilibriumTarget( DomainPartition & domain, // cannot be const...
                                                             std::map< string, localIndex > const & equilNameToEquilId,
                                                             arrayView1d< real64 > const & maxElevation,
                                                             arrayView1d< real64 > const & minElevation ) const
{
  array1d< real64 > localMaxElevation( equilNameToEquilId.size() );
  array1d< real64 > localMinElevation( equilNameToEquilId.size() );
  localMaxElevation.setValues< parallelHostPolicy >( -1e99 );
  localMinElevation.setValues< parallelHostPolicy >( 1e99 );

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply< ElementSubRegionBase,
                   EquilibriumInitialCondition >( 0.0,
                                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                                  EquilibriumInitialCondition::catalogName(),
                                                  [&] ( EquilibriumInitialCondition const & fs,
                                                        string const &,
                                                        SortedArrayView< localIndex const > const & targetSet,
                                                        ElementSubRegionBase & subRegion,
                                                        string const & )
  {
    RAJA::ReduceMax< parallelDeviceReduce, real64 > targetSetMaxElevation( -1e99 );
    RAJA::ReduceMin< parallelDeviceReduce, real64 > targetSetMinElevation( 1e99 );

    arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();

    // TODO: move to FlowSolverBaseKernels to make this function "protected"
    forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const k = targetSet[i];
      targetSetMaxElevation.max( elemCenter[k][2] );
      targetSetMinElevation.min( elemCenter[k][2] );
    } );

    localIndex const equilIndex = equilNameToEquilId.at( fs.getName() );
    localMaxElevation[equilIndex] = LvArray::math::max( targetSetMaxElevation.get(), localMaxElevation[equilIndex] );
    localMinElevation[equilIndex] = LvArray::math::min( targetSetMinElevation.get(), localMinElevation[equilIndex] );

  } );

  MpiWrapper::allReduce( localMaxElevation.toView(),
                         maxElevation,
                         MpiWrapper::Reduction::Max,
                         MPI_COMM_GEOS );
  MpiWrapper::allReduce( localMinElevation.toView(),
                         minElevation,
                         MpiWrapper::Reduction::Min,
                         MPI_COMM_GEOS );
}

void FlowSolverBase::computeSourceFluxSizeScalingFactor( real64 const & time,
                                                         real64 const & dt,
                                                         DomainPartition & domain, // cannot be const...
                                                         std::map< string, localIndex > const & bcNameToBcId,
                                                         arrayView1d< globalIndex > const & bcAllSetsSize ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    fsManager.apply< ElementSubRegionBase,
                     SourceFluxBoundaryCondition >( time + dt,
                                                    mesh,
                                                    SourceFluxBoundaryCondition::catalogName(),
                                                    [&]( SourceFluxBoundaryCondition const & fs,
                                                         string const &,
                                                         SortedArrayView< localIndex const > const & targetSet,
                                                         ElementSubRegionBase & subRegion,
                                                         string const & )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      // TODO: move to FlowSolverBaseKernels to make this function "protected"
      // loop over all the elements of this target set
      RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > localSetSize( 0 );
      forAll< parallelDevicePolicy<> >( targetSet.size(),
                                        [targetSet, ghostRank, localSetSize] GEOS_HOST_DEVICE ( localIndex const k )
      {
        localIndex const ei = targetSet[k];
        if( ghostRank[ei] < 0 )
        {
          localSetSize += 1;
        }
      } );

      // increment the set size for this source flux boundary conditions
      bcAllSetsSize[bcNameToBcId.at( fs.getName())] += localSetSize.get();
    } );
  } );

  // synchronize the set size over all the MPI ranks
  MpiWrapper::allReduce( bcAllSetsSize,
                         bcAllSetsSize,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );
}

void FlowSolverBase::saveAquiferConvergedState( real64 const & time,
                                                real64 const & dt,
                                                DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  // This step requires three passes:
  //    - First we count the number of individual aquifers
  //    - Second we loop over all the stencil entries to compute the sum of aquifer influxes
  //    - Third we loop over the aquifers to save the sums of each individual aquifer

  // Step 1: count individual aquifers

  std::map< string, localIndex > aquiferNameToAquiferId;
  localIndex aquiferCounter = 0;

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
  {
    aquiferNameToAquiferId[bc.getName()] = aquiferCounter;
    aquiferCounter++;
  } );

  // Step 2: sum the aquifer fluxes for each individual aquifer

  array1d< real64 > globalSumFluxes( aquiferNameToAquiferId.size() );
  array1d< real64 > localSumFluxes( aquiferNameToAquiferId.size() );

  fsManager.apply< FaceManager, AquiferBoundaryCondition >( time + dt,
                                                            mesh,
                                                            AquiferBoundaryCondition::catalogName(),
                                                            [&] ( AquiferBoundaryCondition const & bc,
                                                                  string const & setName,
                                                                  SortedArrayView< localIndex const > const &,
                                                                  FaceManager &,
                                                                  string const & )
  {
    BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
    if( stencil.size() == 0 )
    {
      return;
    }

    AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > pressure =
      elemManager.constructFieldAccessor< fields::flow::pressure >();
    pressure.setName( getName() + "/accessors/" + fields::flow::pressure::key() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > pressure_n =
      elemManager.constructFieldAccessor< fields::flow::pressure_n >();
    pressure_n.setName( getName() + "/accessors/" + fields::flow::pressure_n::key() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > gravCoef =
      elemManager.constructFieldAccessor< fields::flow::gravityCoefficient >();
    gravCoef.setName( getName() + "/accessors/" + fields::flow::gravityCoefficient::key() );

    real64 const targetSetSumFluxes = sumAquiferFluxes( stencil,
                                                        aquiferBCWrapper,
                                                        pressure.toNestedViewConst(),
                                                        pressure_n.toNestedViewConst(),
                                                        gravCoef.toNestedViewConst(),
                                                        time,
                                                        dt );

    localIndex const aquiferIndex = aquiferNameToAquiferId.at( bc.getName() );
    localSumFluxes[aquiferIndex] += targetSetSumFluxes;
  } );

  MpiWrapper::allReduce( localSumFluxes,
                         globalSumFluxes,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  // Step 3: we are ready to save the summed fluxes for each individual aquifer

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition & bc )
  {
    localIndex const aquiferIndex = aquiferNameToAquiferId.at( bc.getName() );

    GEOS_LOG_LEVEL_RANK_0_ON_GROUP( logInfo::BoundaryCondition,
                                    GEOS_FMT( "{} {}: at time {} s, the boundary condition produces a volume of {} m3.",
                                              bc.getCatalogName(), bc.getName(),
                                              time + dt, dt * globalSumFluxes[aquiferIndex] ),
                                    bc );
    bc.saveConvergedState( dt * globalSumFluxes[aquiferIndex] );
  } );
}

/**
 * @brief Function to sum the aquiferBC fluxes (as later save them) at the end of the time step
 * This function is applicable for both single-phase and multiphase flow
 */
real64
FlowSolverBase::sumAquiferFluxes( BoundaryStencil const & stencil,
                                  AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & presOld,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  real64 const & timeAtBeginningOfStep,
                                  real64 const & dt )
{
  using Order = BoundaryStencil::Order;

  BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
  BoundaryStencil::WeightContainerViewConstType const & weight = stencil.getWeights();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > targetSetSumFluxes( 0.0 );

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const er  = seri( iconn, Order::ELEM );
    localIndex const esr = sesri( iconn, Order::ELEM );
    localIndex const ei  = sefi( iconn, Order::ELEM );
    real64 const areaFraction = weight( iconn, Order::ELEM );

    // compute the aquifer influx rate using the pressure influence function and the aquifer props
    real64 dAquiferVolFlux_dPres = 0.0;
    real64 const aquiferVolFlux = aquiferBCWrapper.compute( timeAtBeginningOfStep,
                                                            dt,
                                                            pres[er][esr][ei],
                                                            presOld[er][esr][ei],
                                                            gravCoef[er][esr][ei],
                                                            areaFraction,
                                                            dAquiferVolFlux_dPres );
    targetSetSumFluxes += aquiferVolFlux;
  } );
  return targetSetSumFluxes.get();
}

void FlowSolverBase::prepareStencilWeights( DomainPartition & domain ) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( getDiscretizationName() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > hydraulicAperture =
      m_isLaggingFractureStencilWeightsUpdate ?
      mesh.getElemManager().constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::aperture0::key() ) :
      mesh.getElemManager().constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::hydraulicAperture::key() );

    fluxApprox.forStencils< SurfaceElementStencil, FaceElementToCellStencil, EmbeddedSurfaceToCellStencil >( mesh, [&]( auto & stencil )
    {
      using STENCILWRAPPER_TYPE =  typename TYPEOFREF( stencil ) ::KernelWrapper;

      STENCILWRAPPER_TYPE stencilWrapper = stencil.createKernelWrapper();

      flowSolverBaseKernels::stencilWeightsUpdateKernel< STENCILWRAPPER_TYPE >::prepareStencilWeights( stencilWrapper, hydraulicAperture.toNestedViewConst() );
    } );
  } );
}

void FlowSolverBase::updateStencilWeights( DomainPartition & domain ) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( getDiscretizationName() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > hydraulicAperture =
      mesh.getElemManager().constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::hydraulicAperture::key() );

    fluxApprox.forStencils< SurfaceElementStencil, FaceElementToCellStencil, EmbeddedSurfaceToCellStencil >( mesh, [&]( auto & stencil )
    {
      using STENCILWRAPPER_TYPE =  typename TYPEOFREF( stencil ) ::KernelWrapper;

      STENCILWRAPPER_TYPE stencilWrapper = stencil.createKernelWrapper();

      flowSolverBaseKernels::stencilWeightsUpdateKernel< STENCILWRAPPER_TYPE >::updateStencilWeights( stencilWrapper, hydraulicAperture.toNestedViewConst() );
    } );
  } );
}

bool FlowSolverBase::checkSequentialSolutionIncrements( DomainPartition & GEOS_UNUSED_PARAM( domain ) ) const
{

  GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                         GEOS_FMT( "    {}: Max pressure change during outer iteration: {} Pa",
                                   getName(), GEOS_FMT( "{:.{}f}", m_sequentialPresChange, 3 ) ) );

  if( m_isThermal )
  {
    GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                           GEOS_FMT( "    {}: Max temperature change during outer iteration: {} K",
                                     getName(), GEOS_FMT( "{:.{}f}", m_sequentialTempChange, 3 ) ) );
  }

  return (m_sequentialPresChange < m_maxSequentialPresChange) && (m_sequentialTempChange < m_maxSequentialTempChange);
}

string FlowSolverBase::BCMessage::generateMessage( string_view baseMessage,
                                                   string_view fieldName, string_view setName )
{
  return GEOS_FMT( "{}\nCheck if you have added or applied the appropriate fields to "
                   "the FieldSpecification component with fieldName=\"{}\" "
                   "and setNames=\"{}\"\n", baseMessage, fieldName, setName );
}

string FlowSolverBase::BCMessage::pressureConflict( string_view regionName, string_view subRegionName,
                                                    string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Conflicting pressure boundary conditions on set {}/{}/{}",
                                    regionName, subRegionName, setName ), fieldName, setName );
}

string FlowSolverBase::BCMessage::temperatureConflict( string_view regionName, string_view subRegionName,
                                                       string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Conflicting temperature boundary conditions on set {}/{}/{}",
                                    regionName, subRegionName, setName ), fieldName, setName );
}

string FlowSolverBase::BCMessage::missingPressure( string_view regionName, string_view subRegionName,
                                                   string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Pressure boundary condition not prescribed on set {}/{}/{}",
                                    regionName, subRegionName, setName ), fieldName, setName );
}

string FlowSolverBase::BCMessage::missingTemperature( string_view regionName, string_view subRegionName,
                                                      string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Temperature boundary condition not prescribed on set {}/{}/{}",
                                    regionName, subRegionName, setName ), fieldName, setName );
}

string FlowSolverBase::BCMessage::conflictingComposition( int comp, string_view componentName,
                                                          string_view regionName, string_view subRegionName,
                                                          string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Conflicting {} composition (no.{}) for boundary conditions on set {}/{}/{}",
                                    componentName, comp, regionName, subRegionName, setName ),
                          fieldName, setName );
}

string FlowSolverBase::BCMessage::invalidComponentIndex( int comp,
                                                         string_view fsName,
                                                         string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Invalid component index no.{} in boundary condition {}",
                                    comp, fsName ), fieldName, fsName );
}

string FlowSolverBase::BCMessage::notAppliedOnRegion( int componentIndex, string_view componentName,
                                                      string_view regionName, string_view subRegionName,
                                                      string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Boundary condition not applied to {} component (no.{})"
                                    "on region {}/{}/{}\n",
                                    componentName, componentIndex, regionName, subRegionName, setName ),
                          fieldName, setName );
}

} // namespace geos
