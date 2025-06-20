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
 * @file CompositionalMultiphaseFVM.cpp
 */

#include "CompositionalMultiphaseFVM.hpp"

#include "common/MpiWrapper.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "dataRepository/Group.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/LogLevelsInfo.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ResidualNormKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalResidualNormKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalSolutionScalingKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalSolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/FluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/DiffusionDispersionFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalDiffusionDispersionFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/StabilizedFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/DissipationFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/DirichletFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalDirichletFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/PhaseMobilityKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalPhaseMobilityKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/AquiferBCKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/CFLKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/zFormulation/FluxComputeZFormulationKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/zFormulation/AccumulationZFormulationKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/zFormulation/PhaseMobilityZFormulationKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/zFormulation/SolutionScalingZFormulationKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/zFormulation/DirichletFluxComputeZFormulationKernel.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/MultiphasePoromechanicsConformingFractures.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace compositionalMultiphaseUtilities; // for ScalingType

CompositionalMultiphaseFVM::CompositionalMultiphaseFVM( const string & name,
                                                        Group * const parent )
  :
  CompositionalMultiphaseBase( name, parent )
{
  registerWrapper( viewKeyStruct::useDBCString(), &m_dbcParams.useDBC ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Enable Dissipation-based continuation flux" );

  registerWrapper( viewKeyStruct::omegaDBCString(), &m_dbcParams.omega ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which DBC flux is multiplied" );

  registerWrapper( viewKeyStruct::continuationDBCString(), &m_dbcParams.continuation ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for enabling continuation parameter" );

  registerWrapper( viewKeyStruct::miscibleDBCString(), &m_dbcParams.miscible ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for enabling DBC formulation with/without miscibility" );

  registerWrapper( viewKeyStruct::kappaminDBCString(), &m_dbcParams.kappamin ).
    setApplyDefaultValue( 1e-20 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor that controls how much dissipation is kept in the system when continuation is used" );

  registerWrapper( viewKeyStruct::contMultiplierDBCString(), &m_dbcParams.contMultiplier ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which continuation parameter is changed every newton when DBC is used" );

  registerWrapper( viewKeyStruct::scalingTypeString(), &m_scalingType ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( ScalingType::Global ).
    setDescription( "Solution scaling type."
                    "Valid options:\n* " + EnumStrings< ScalingType >::concat( "\n* " ) );

  this->registerWrapper( viewKeyStruct::gravityDensitySchemeString(), &m_gravityDensityScheme ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( GravityDensityScheme::ArithmeticAverage ).
    setDescription( "Scheme for density treatment in gravity" );

  this->registerWrapper( viewKeyStruct::targetFlowCFLString(), &m_targetFlowCFL ).
    setApplyDefaultValue( -1. ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target CFL condition `CFL condition <http://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition>`_"
                    " when computing the next timestep." );

  addLogLevel< logInfo::Convergence >();
  addLogLevel< logInfo::Solution >();
  addLogLevel< logInfo::TimeStep >();
}

void CompositionalMultiphaseFVM::postInputInitialization()
{
  CompositionalMultiphaseBase::postInputInitialization();

  if( m_scalingType == ScalingType::Local &&
      m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None )
  {
    GEOS_ERROR( GEOS_FMT( "{}: line search is not supported for {} = {}",
                          getName(), viewKeyStruct::scalingTypeString(),
                          EnumStrings< ScalingType >::toString( ScalingType::Local )) );
  }

  if( m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
  {
    string const formulationName = EnumStrings< CompositionalMultiphaseFormulationType >::toString( CompositionalMultiphaseFormulationType::OverallComposition );

    if( m_dbcParams.useDBC ) // z_c formulation is not compatible with DBC
    {
      GEOS_ERROR( GEOS_FMT( "{}: '{}' is not compatible with {}",
                            getDataContext(), formulationName, viewKeyStruct::useDBCString() ) );
    }

    DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
    auto const & upwindingParams = fluxApprox.upwindingParams();
    if( upwindingParams.upwindingScheme == UpwindingScheme::C1PPU ||
        upwindingParams.upwindingScheme == UpwindingScheme::IHU )
    {
      GEOS_ERROR( GEOS_FMT( "{}: {} is not available for {}",
                            getDataContext(),
                            EnumStrings< UpwindingScheme >::toString( upwindingParams.upwindingScheme ),
                            formulationName ) );
    }
  }
}

void CompositionalMultiphaseFVM::registerDataOnMesh( Group & meshBodies )
{
  using namespace fields::flow;

  CompositionalMultiphaseBase::registerDataOnMesh( meshBodies );

  if( m_targetFlowCFL > 0 )
  {
    registerDataForCFL( meshBodies );
  }
}

void CompositionalMultiphaseFVM::registerDataForCFL( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< fields::flow::phaseOutflux >( getName()).reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< fields::flow::componentOutflux >( getName()).reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< fields::flow::phaseCFLNumber >( getName());
      subRegion.registerField< fields::flow::componentCFLNumber >( getName());
    } );
  } );
}

void CompositionalMultiphaseFVM::initializePreSubGroups()
{
  CompositionalMultiphaseBase::initializePreSubGroups();

  m_linearSolverParameters.get().mgr.strategy = m_isThermal
                                                ? LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseFVM
                                                : LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM;

  checkDiscretizationName();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  GEOS_ERROR_IF( fluxApprox.upwindingParams().upwindingScheme == UpwindingScheme::HU2PH && m_numPhases != 2,
                 GEOS_FMT( "{}: upwinding scheme {} only supports 2-phase flow",
                           getName(), EnumStrings< UpwindingScheme >::toString( UpwindingScheme::HU2PH )));
}

void CompositionalMultiphaseFVM::setupDofs( DomainPartition const & domain,
                                            DofManager & dofManager ) const
{
  // add a field for the cell-centered degrees of freedom
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  // this call with instruct GEOS to reorder the dof numbers
  dofManager.setLocalReorderingType( viewKeyStruct::elemDofFieldString(),
                                     DofManager::LocalReorderingType::ReverseCutHillMcKee );

  // for the volume balance equation, disable global coupling
  // this equation is purely local (not coupled to neighbors or other physics)
  dofManager.disableGlobalCouplingForEquation( viewKeyStruct::elemDofFieldString(),
                                               m_numComponents );


  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}


void CompositionalMultiphaseFVM::assembleFluxTerms( real64 const dt,
                                                    DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  using namespace isothermalCompositionalMultiphaseFVMKernels;

  BitFlags< KernelFlags > kernelFlags;
  if( m_hasCapPressure )
    kernelFlags.set( KernelFlags::CapPressure );
  if( m_hasDiffusion )
    kernelFlags.set( KernelFlags::Diffusion );
  if( m_hasDispersion )
    kernelFlags.set( KernelFlags::Dispersion );
  if( m_useTotalMassEquation )
    kernelFlags.set( KernelFlags::TotalMassEquation );
  if( m_gravityDensityScheme == GravityDensityScheme::PhasePresence )
    kernelFlags.set( KernelFlags::CheckPhasePresenceInGravity );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    auto const & upwindingParams = fluxApprox.upwindingParams();
    if( upwindingParams.upwindingScheme == UpwindingScheme::C1PPU &&
        isothermalCompositionalMultiphaseFVMKernelUtilities::epsC1PPU > 0 )
      kernelFlags.set( KernelFlags::C1PPU );
    else if( upwindingParams.upwindingScheme == UpwindingScheme::IHU )
      kernelFlags.set( KernelFlags::IHU );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      // Convective flux

      if( m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
      {
        // isothermal only for now
        isothermalCompositionalMultiphaseFVMKernels::
          FluxComputeZFormulationKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     elemDofKey,
                                                     kernelFlags,
                                                     getName(),
                                                     mesh.getElemManager(),
                                                     stencilWrapper,
                                                     dt,
                                                     localMatrix.toViewConstSizes(),
                                                     localRhs.toView() );
      }
      else
      {
        if( m_isThermal )
        {
          thermalCompositionalMultiphaseFVMKernels::
            FluxComputeKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       dofManager.rankOffset(),
                                                       elemDofKey,
                                                       kernelFlags,
                                                       getName(),
                                                       mesh.getElemManager(),
                                                       stencilWrapper,
                                                       dt,
                                                       localMatrix.toViewConstSizes(),
                                                       localRhs.toView() );
        }
        else
        {
          if( m_dbcParams.useDBC )
          {
            dissipationCompositionalMultiphaseFVMKernels::
              FluxComputeKernelFactory::
              createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                         m_numPhases,
                                                         dofManager.rankOffset(),
                                                         elemDofKey,
                                                         kernelFlags,
                                                         getName(),
                                                         mesh.getElemManager(),
                                                         stencilWrapper,
                                                         dt,
                                                         localMatrix.toViewConstSizes(),
                                                         localRhs.toView(),
                                                         m_dbcParams.omega,
                                                         getNonlinearSolverParameters().m_numNewtonIterations,
                                                         m_dbcParams.continuation,
                                                         m_dbcParams.miscible,
                                                         m_dbcParams.kappamin,
                                                         m_dbcParams.contMultiplier );
          }
          else
          {
            isothermalCompositionalMultiphaseFVMKernels::
              FluxComputeKernelFactory::
              createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                         m_numPhases,
                                                         dofManager.rankOffset(),
                                                         elemDofKey,
                                                         kernelFlags,
                                                         getName(),
                                                         mesh.getElemManager(),
                                                         stencilWrapper,
                                                         dt,
                                                         localMatrix.toViewConstSizes(),
                                                         localRhs.toView() );
          }
        }

        // Diffusive and dispersive flux

        if( m_hasDiffusion || m_hasDispersion )
        {

          if( m_isThermal )
          {
            thermalCompositionalMultiphaseFVMKernels::
              DiffusionDispersionFluxComputeKernelFactory::
              createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                         m_numPhases,
                                                         dofManager.rankOffset(),
                                                         elemDofKey,
                                                         kernelFlags,
                                                         getName(),
                                                         mesh.getElemManager(),
                                                         stencilWrapper,
                                                         dt,
                                                         localMatrix.toViewConstSizes(),
                                                         localRhs.toView() );
          }
          else
          {
            isothermalCompositionalMultiphaseFVMKernels::
              DiffusionDispersionFluxComputeKernelFactory::
              createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                         m_numPhases,
                                                         dofManager.rankOffset(),
                                                         elemDofKey,
                                                         kernelFlags,
                                                         getName(),
                                                         mesh.getElemManager(),
                                                         stencilWrapper,
                                                         dt,
                                                         localMatrix.toViewConstSizes(),
                                                         localRhs.toView() );
          }
        }
      }

    } );
  } );
}

void CompositionalMultiphaseFVM::assembleStabilizedFluxTerms( real64 const dt,
                                                              DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  using namespace isothermalCompositionalMultiphaseFVMKernels;

  BitFlags< KernelFlags > kernelFlags;
  if( m_hasCapPressure )
    kernelFlags.set( KernelFlags::CapPressure );
  if( m_useTotalMassEquation )
    kernelFlags.set( KernelFlags::TotalMassEquation );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      // Thermal implementation not supported yet

      stabilizedCompositionalMultiphaseFVMKernels::
        FluxComputeKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   m_numPhases,
                                                   dofManager.rankOffset(),
                                                   elemDofKey,
                                                   kernelFlags,
                                                   getName(),
                                                   mesh.getElemManager(),
                                                   stencilWrapper,
                                                   dt,
                                                   localMatrix.toViewConstSizes(),
                                                   localRhs.toView() );

    } );
  } );
}

real64 CompositionalMultiphaseFVM::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                          real64 const & GEOS_UNUSED_PARAM( dt ),
                                                          DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  integer constexpr numNorm = 3; // mass/volume balance and energy balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );

  physicsSolverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[numNorm]{};
      real64 subRegionResidualNormalizer[numNorm]{};

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      // step 1: compute the norm in the subRegion

      if( m_isThermal )
      {
        string const & solidInternalEnergyName = subRegion.getReference< string >( viewKeyStruct::solidInternalEnergyNamesString() );
        SolidInternalEnergy const & solidInternalEnergy = getConstitutiveModel< SolidInternalEnergy >( subRegion, solidInternalEnergyName );

        thermalCompositionalMultiphaseBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     numFluidComponents(),
                                                     numFluidPhases(),
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     solidInternalEnergy,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionResidualNorm,
                                                     subRegionResidualNormalizer );
      }
      else
      {
        real64 subRegionFlowResidualNorm[2]{};
        real64 subRegionFlowResidualNormalizer[2]{};
        isothermalCompositionalMultiphaseBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     numFluidComponents(),
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionFlowResidualNorm,
                                                     subRegionFlowResidualNormalizer );
        // mass
        subRegionResidualNorm[0] = subRegionFlowResidualNorm[0];
        subRegionResidualNormalizer[0] = subRegionFlowResidualNormalizer[0];
        // volume
        subRegionResidualNorm[1] = subRegionFlowResidualNorm[1];
        subRegionResidualNormalizer[1] = subRegionFlowResidualNormalizer[1];
      }

      // step 2: first reduction across meshBodies/regions/subRegions

      if( normType == physicsSolverBaseKernels::NormType::Linf )
      {
        physicsSolverBaseKernels::LinfResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, localResidualNorm );
      }
      else
      {
        physicsSolverBaseKernels::L2ResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, subRegionResidualNormalizer, localResidualNorm, localResidualNormalizer );
      }
    } );
  } );

  // step 3: second reduction across MPI ranks

  real64 residualNorm = 0.0;
  array1d< real64 > globalResidualNorm;
  globalResidualNorm.resize( numNorm );
  if( m_isThermal )
  {
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      physicsSolverBaseKernels::LinfResidualNormHelper::
        computeGlobalNorm( localResidualNorm, globalResidualNorm );
    }
    else
    {
      physicsSolverBaseKernels::L2ResidualNormHelper::
        computeGlobalNorm( localResidualNorm, localResidualNormalizer, globalResidualNorm );
    }
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1]  + globalResidualNorm[2] * globalResidualNorm[2] );

    GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::Convergence,
                               GEOS_FMT( "        ( Rmass Rvol ) = ( {:4.2e} {:4.2e} )        ( Renergy ) = ( {:4.2e} )",
                                         globalResidualNorm[0], globalResidualNorm[1], globalResidualNorm[2] ));
  }
  else
  {
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      physicsSolverBaseKernels::LinfResidualNormHelper::
        computeGlobalNorm( localResidualNorm, globalResidualNorm );
    }
    else
    {
      physicsSolverBaseKernels::L2ResidualNormHelper::
        computeGlobalNorm( localResidualNorm, localResidualNormalizer, globalResidualNorm );
    }
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );

    GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::Convergence, GEOS_FMT( "        ( Rmass Rvol ) = ( {:4.2e} {:4.2e} )",
                                                               globalResidualNorm[0], globalResidualNorm[1] ) );
  }

  return residualNorm;
}

real64 CompositionalMultiphaseFVM::scalingForSystemSolution( DomainPartition & domain,
                                                             DofManager const & dofManager,
                                                             arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  if( m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
  {
    return scalingForSystemSolutionZFormulation( domain, dofManager, localSolution );
  }
  else
  {
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    real64 scalingFactor = 1.0;
    real64 minPresScalingFactor = 1.0, minCompDensScalingFactor = 1.0, minTempScalingFactor = 1.0;

    stdVector< MpiWrapper::PairType< real64, globalIndex > > regionDeltaPresMaxLoc;
    stdVector< MpiWrapper::PairType< real64, globalIndex > > regionDeltaCompDensMaxLoc;
    stdVector< MpiWrapper::PairType< real64, globalIndex > > regionDeltaTempMaxLoc;

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        arrayView1d< globalIndex const > const localToGlobalMap = subRegion.localToGlobalMap();
        arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
        arrayView1d< real64 const > const temperature = subRegion.getField< fields::flow::temperature >();
        arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::flow::globalCompDensity >();
        arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
        arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::flow::temperatureScalingFactor >();
        arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::flow::globalCompDensityScalingFactor >();

        const integer temperatureOffset = m_numComponents+1;

        auto const subRegionData = m_isThermal ?
                                   thermalCompositionalMultiphaseBaseKernels::
                                     SolutionScalingKernelFactory::
                                     createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                                                m_maxAbsolutePresChange,
                                                                                m_maxRelativeTempChange,
                                                                                m_maxCompFracChange,
                                                                                m_maxRelativeCompDensChange,
                                                                                pressure,
                                                                                temperature,
                                                                                compDens,
                                                                                pressureScalingFactor,
                                                                                compDensScalingFactor,
                                                                                temperatureScalingFactor,
                                                                                dofManager.rankOffset(),
                                                                                m_numComponents,
                                                                                dofKey,
                                                                                subRegion,
                                                                                localSolution,
                                                                                temperatureOffset ):
                                   isothermalCompositionalMultiphaseBaseKernels::
                                     SolutionScalingKernelFactory::
                                     createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                                                m_maxAbsolutePresChange,
                                                                                m_maxCompFracChange,
                                                                                m_maxRelativeCompDensChange,
                                                                                pressure,
                                                                                compDens,
                                                                                pressureScalingFactor,
                                                                                compDensScalingFactor,
                                                                                dofManager.rankOffset(),
                                                                                m_numComponents,
                                                                                dofKey,
                                                                                subRegion,
                                                                                localSolution );

        if( subRegion.size() > 0 || subRegion.size() !=  subRegion.getNumberOfGhosts() )
        {
          if( m_scalingType == ScalingType::Global )
          {
            scalingFactor = std::min( scalingFactor, subRegionData.localMinVal );
          }

          regionDeltaPresMaxLoc.push_back( { subRegionData.localMaxDeltaPres,
                                             subRegionData.localMaxDeltaPresLoc >= 0 ? localToGlobalMap[subRegionData.localMaxDeltaPresLoc] : -1 } );
          minPresScalingFactor = std::min( minPresScalingFactor, subRegionData.localMinPresScalingFactor );

          regionDeltaCompDensMaxLoc.push_back( { subRegionData.localMaxDeltaCompDens,
                                                 subRegionData.localMaxDeltaCompDensLoc >= 0 ? localToGlobalMap[subRegionData.localMaxDeltaCompDensLoc] : -1 } );
          minCompDensScalingFactor = std::min( minCompDensScalingFactor, subRegionData.localMinCompDensScalingFactor );

          if( m_isThermal )
          {
            regionDeltaTempMaxLoc.push_back( { subRegionData.localMaxDeltaTemp,
                                               subRegionData.localMaxDeltaTempLoc >= 0 ? localToGlobalMap[subRegionData.localMaxDeltaTempLoc] : -1 } );
            minTempScalingFactor = std::min( minTempScalingFactor, subRegionData.localMinTempScalingFactor );
          }
        }
      } );
    } );

    auto globalDeltaPresMax = MpiWrapper::max< real64, globalIndex >( regionDeltaPresMaxLoc );
    auto globalDeltaCompDensMax = MpiWrapper::max< real64, globalIndex >( regionDeltaCompDensMaxLoc );

    scalingFactor = MpiWrapper::min( scalingFactor );
    minPresScalingFactor = MpiWrapper::min( minPresScalingFactor );
    minCompDensScalingFactor = MpiWrapper::min( minCompDensScalingFactor );

    string const massUnit = m_useMass ? "kg/m3" : "mol/m3";
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Max pressure change = {:.3f} Pa (before scaling) at cell {}",
                                     getName(), globalDeltaPresMax.first, globalDeltaPresMax.second ) );

    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Max component density change = {:.3f} {} (before scaling) at cell {}",
                                     getName(), globalDeltaCompDensMax.first, massUnit, globalDeltaCompDensMax.second ) );

    if( m_isThermal )
    {
      auto globalMaxDeltaTemp = MpiWrapper::max< real64, globalIndex >( regionDeltaTempMaxLoc );

      minTempScalingFactor = MpiWrapper::min( minTempScalingFactor );
      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                             GEOS_FMT( "        {}: Max temperature change = {:.3f} K (before scaling) at cell maxRegionDeltaTempLoc {}",
                                       getName(), globalMaxDeltaTemp.first, globalMaxDeltaTemp.second ) );
    }

    if( m_scalingType == ScalingType::Local )
    {
      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min pressure scaling factor = {}", getName(), minPresScalingFactor ) );
      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min component density scaling factor = {}", getName(), minCompDensScalingFactor ) );
      if( m_isThermal )
      {
        GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min temperature scaling factor = {}", getName(), minTempScalingFactor ) );
      }
    }

    return LvArray::math::max( scalingFactor, m_minScalingFactor );
  }
}

real64 CompositionalMultiphaseFVM::scalingForSystemSolutionZFormulation( DomainPartition & domain,
                                                                         DofManager const & dofManager,
                                                                         arrayView1d< real64 const > const & localSolution )
{
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  real64 scalingFactor = 1.0;
  real64 maxDeltaPres = 0.0, maxDeltaCompFrac = 0.0, maxDeltaTemp = 0.0;
  real64 minPresScalingFactor = 1.0, minCompFracScalingFactor = 1.0, minTempScalingFactor = 1.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
      //arrayView1d< real64 const > const temperature = subRegion.getField< fields::flow::temperature >();
      arrayView2d< real64 const, compflow::USD_COMP > const compFrac = subRegion.getField< fields::flow::globalCompFraction >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      //arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::flow::temperatureScalingFactor >();
      arrayView1d< real64 > compFracScalingFactor = subRegion.getField< fields::flow::globalCompFractionScalingFactor >();

      auto const subRegionData = isothermalCompositionalMultiphaseBaseKernels::
                                   SolutionScalingZFormulationKernelFactory::
                                   createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                                              m_maxAbsolutePresChange,
                                                                              m_maxCompFracChange,
                                                                              pressure,
                                                                              compFrac,
                                                                              pressureScalingFactor,
                                                                              compFracScalingFactor,
                                                                              dofManager.rankOffset(),
                                                                              m_numComponents,
                                                                              dofKey,
                                                                              subRegion,
                                                                              localSolution );

      if( m_scalingType == ScalingType::Global )
      {
        scalingFactor = std::min( scalingFactor, subRegionData.localMinVal );
      }
      maxDeltaPres  = std::max( maxDeltaPres, subRegionData.localMaxDeltaPres );
      maxDeltaCompFrac = std::max( maxDeltaCompFrac, subRegionData.localMaxDeltaCompFrac );
      maxDeltaTemp = std::max( maxDeltaTemp, subRegionData.localMaxDeltaTemp );
      minPresScalingFactor = std::min( minPresScalingFactor, subRegionData.localMinPresScalingFactor );
      minCompFracScalingFactor = std::min( minCompFracScalingFactor, subRegionData.localMinCompFracScalingFactor );
      minTempScalingFactor = std::min( minTempScalingFactor, subRegionData.localMinTempScalingFactor );

    } );
  } );

  scalingFactor = MpiWrapper::min( scalingFactor );
  maxDeltaPres  = MpiWrapper::max( maxDeltaPres );
  maxDeltaCompFrac = MpiWrapper::max( maxDeltaCompFrac );
  minPresScalingFactor = MpiWrapper::min( minPresScalingFactor );
  minCompFracScalingFactor = MpiWrapper::min( minCompFracScalingFactor );

  GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Max pressure change = {} Pa (before scaling)",
                                                      getName(), GEOS_FMT( "{:.{}f}", maxDeltaPres, 3 ) ) );
  GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Max component fraction change = {} (before scaling)",
                                                      getName(), GEOS_FMT( "{:.{}f}", maxDeltaCompFrac, 3 ) ) );

  if( m_isThermal )
  {
    maxDeltaTemp = MpiWrapper::max( maxDeltaTemp );
    minTempScalingFactor = MpiWrapper::min( minTempScalingFactor );
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Max temperature change = {} K (before scaling)",
                                                        getName(), GEOS_FMT( "{:.{}f}", maxDeltaTemp, 3 ) ) );
  }

  if( m_scalingType == ScalingType::Local )
  {
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min pressure scaling factor = {}", getName(), minPresScalingFactor ) );
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min component fraction scaling factor = {}", getName(), minCompFracScalingFactor ) );
    if( m_isThermal )
    {
      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                             GEOS_FMT( "        {}: Min temperature scaling factor = {}", getName(), minTempScalingFactor ) );
    }
  }

  return LvArray::math::max( scalingFactor, m_minScalingFactor );
}

bool CompositionalMultiphaseFVM::checkSystemSolution( DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution,
                                                      real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  if( m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
  {
    // TO DO: Implement the solution check for Z Formulation
    return true;
  }
  else
  {
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    integer localCheck = 1;
    real64 minPres = 0.0, minDens = 0.0, minTotalDens = 0.0;
    integer numNegPres = 0, numNegDens = 0, numNegTotalDens = 0;

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        arrayView1d< real64 const > const pressure =
          subRegion.getField< fields::flow::pressure >();
        arrayView1d< real64 const > const temperature =
          subRegion.getField< fields::flow::temperature >();
        arrayView2d< real64 const, compflow::USD_COMP > const compDens =
          subRegion.getField< fields::flow::globalCompDensity >();
        arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
        arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::flow::temperatureScalingFactor >();
        arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::flow::globalCompDensityScalingFactor >();
        // check that pressure and component densities are non-negative
        // for thermal, check that temperature is above 273.15 K
        const integer temperatureOffset = m_numComponents+1;
        auto const subRegionData =
          m_isThermal
    ? thermalCompositionalMultiphaseBaseKernels::
            SolutionCheckKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                       m_allowNegativePressure,
                                                       m_scalingType,
                                                       scalingFactor,
                                                       pressure,
                                                       temperature,
                                                       compDens,
                                                       pressureScalingFactor,
                                                       temperatureScalingFactor,
                                                       compDensScalingFactor,
                                                       dofManager.rankOffset(),
                                                       m_numComponents,
                                                       dofKey,
                                                       subRegion,
                                                       localSolution,
                                                       temperatureOffset )
    : isothermalCompositionalMultiphaseBaseKernels::
            SolutionCheckKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                       m_allowNegativePressure,
                                                       m_scalingType,
                                                       scalingFactor,
                                                       pressure,
                                                       compDens,
                                                       pressureScalingFactor,
                                                       compDensScalingFactor,
                                                       dofManager.rankOffset(),
                                                       m_numComponents,
                                                       dofKey,
                                                       subRegion,
                                                       localSolution );

        localCheck = std::min( localCheck, subRegionData.localMinVal );

        minPres  = std::min( minPres, subRegionData.localMinPres );
        minDens = std::min( minDens, subRegionData.localMinDens );
        minTotalDens = std::min( minTotalDens, subRegionData.localMinTotalDens );
        numNegPres += subRegionData.localNumNegPressures;
        numNegDens += subRegionData.localNumNegDens;
        numNegTotalDens += subRegionData.localNumNegTotalDens;
      } );
    } );

    minPres  = MpiWrapper::min( minPres );
    minDens = MpiWrapper::min( minDens );
    minTotalDens = MpiWrapper::min( minTotalDens );
    numNegPres = MpiWrapper::sum( numNegPres );
    numNegDens = MpiWrapper::sum( numNegDens );
    numNegTotalDens = MpiWrapper::sum( numNegTotalDens );

    if( numNegPres > 0 )
      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                             GEOS_FMT( "        {}: Number of negative pressure values: {}, minimum value: {} Pa",
                                       getName(), numNegPres, fmt::format( "{:.{}f}", minPres, 3 ) ) );
    string const massUnit = m_useMass ? "kg/m3" : "mol/m3";
    if( numNegDens > 0 )
      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                             GEOS_FMT( "        {}: Number of negative component density values: {}, minimum value: {} {}}",
                                       getName(), numNegDens, fmt::format( "{:.{}f}", minDens, 3 ), massUnit ) );
    if( minTotalDens > 0 )
      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                             GEOS_FMT( "        {}: Number of negative total density values: {}, minimum value: {} {}}",
                                       getName(), minTotalDens, fmt::format( "{:.{}f}", minDens, 3 ), massUnit ) );

    return MpiWrapper::min( localCheck );
  }
}

void CompositionalMultiphaseFVM::applySystemSolution( DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution,
                                                      real64 const scalingFactor,
                                                      real64 const GEOS_UNUSED_PARAM( dt ),
                                                      DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  bool const localScaling = m_scalingType == ScalingType::Local;

  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
  DofManager::CompMask componentMask( m_numDofPerCell, 1, m_numComponents+1 );

  if( localScaling )
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 fields::flow::pressureScalingFactor::key(),
                                 pressureMask );
  }
  else
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor,
                                 pressureMask );
  }

  if( m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
  {
    if( localScaling )
    {
      dofManager.addVectorToField( localSolution,
                                   viewKeyStruct::elemDofFieldString(),
                                   fields::flow::globalCompFraction::key(),
                                   fields::flow::globalCompFractionScalingFactor::key(),
                                   componentMask );
    }
    else
    {
      dofManager.addVectorToField( localSolution,
                                   viewKeyStruct::elemDofFieldString(),
                                   fields::flow::globalCompFraction::key(),
                                   scalingFactor,
                                   componentMask );
    }
  }
  else
  {
    if( localScaling )
    {
      dofManager.addVectorToField( localSolution,
                                   viewKeyStruct::elemDofFieldString(),
                                   fields::flow::globalCompDensity::key(),
                                   fields::flow::globalCompDensityScalingFactor::key(),
                                   componentMask );
    }
    else
    {
      dofManager.addVectorToField( localSolution,
                                   viewKeyStruct::elemDofFieldString(),
                                   fields::flow::globalCompDensity::key(),
                                   scalingFactor,
                                   componentMask );
    }
  }

  if( m_isThermal )
  {
    DofManager::CompMask temperatureMask( m_numDofPerCell, m_numComponents+1, m_numComponents+2 );
    if( localScaling )
    {
      dofManager.addVectorToField( localSolution,
                                   viewKeyStruct::elemDofFieldString(),
                                   fields::flow::temperature::key(),
                                   fields::flow::temperatureScalingFactor::key(),
                                   temperatureMask );
    }
    else
    {
      dofManager.addVectorToField( localSolution,
                                   viewKeyStruct::elemDofFieldString(),
                                   fields::flow::temperature::key(),
                                   scalingFactor,
                                   temperatureMask );
    }
  }

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function

  if( m_allowCompDensChopping )
  {
    if( m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
      chopNegativeCompFractions( domain );
    else
      chopNegativeDensities( domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    stdVector< string > fields{ fields::flow::pressure::key(),
                                m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition ?
                                fields::flow::globalCompFraction::key() : fields::flow::globalCompDensity::key() };
    if( m_isThermal )
    {
      fields.emplace_back( fields::flow::temperature::key() );
    }
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void CompositionalMultiphaseFVM::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  // note that the phase mobility computed here also includes phase density
  string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, fluidName );

  string const & relpermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase const & relperm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relpermName );

  if( m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
  {
    // For now: isothermal only
    isothermalCompositionalMultiphaseFVMKernels::
      PhaseMobilityZFormulationKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 dataGroup,
                                                 fluid,
                                                 relperm );
  }
  else
  {
    if( m_isThermal )
    {
      thermalCompositionalMultiphaseFVMKernels::
        PhaseMobilityKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   m_numPhases,
                                                   dataGroup,
                                                   fluid,
                                                   relperm );
    }
    else
    {
      isothermalCompositionalMultiphaseFVMKernels::
        PhaseMobilityKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   m_numPhases,
                                                   dataGroup,
                                                   fluid,
                                                   relperm );
    }
  }
}

void CompositionalMultiphaseFVM::applyBoundaryConditions( real64 time_n,
                                                          real64 dt,
                                                          DomainPartition & domain,
                                                          DofManager const & dofManager,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  CompositionalMultiphaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
  if( !m_keepVariablesConstantDuringInitStep )
  {
    applyFaceDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
  }
}

bool CompositionalMultiphaseFVM::validateFaceDirichletBC( DomainPartition & domain,
                                                          real64 const time ) const
{
  constexpr integer MAX_NC = MultiFluidBase::MAX_NUM_COMPONENTS;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  bool bcConsistent = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {

    // maps to check consistent application of BC
    // maps: setName (-> numComps)
    map< string, ComponentMask< MAX_NC > > bcPresCompStatusMap; // check that pressure/comp are present/consistent
    set< string > bcTempStatusMap; // check that temperature is present/consistent

    // 1. Check pressure Dirichlet BCs
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    fields::flow::pressure::key(),
                                    [&]( FieldSpecificationBase const &,
                                         string const & setName,
                                         SortedArrayView< localIndex const > const &,
                                         FaceManager &,
                                         string const & )
    {
      // Check whether pressure has already been applied to this set
      if( bcPresCompStatusMap.count( setName ) > 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Conflicting pressure boundary conditions on set {}", setName ) );
      }
      bcPresCompStatusMap[setName].setNumComp( m_numComponents );
    } );

    // 2. Check temperature Dirichlet BCs (we always require a temperature for face-based BCs)
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    fields::flow::temperature::key(),
                                    [&]( FieldSpecificationBase const &,
                                         string const & setName,
                                         SortedArrayView< localIndex const > const &,
                                         FaceManager &,
                                         string const & )
    {
      // 2.1 Check whether temperature has already been applied to this set
      if( bcTempStatusMap.count( setName ) > 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Conflicting temperature boundary conditions on set {}", setName ) );
      }
      bcTempStatusMap.insert( setName );

      // 2.2 Check that there is pressure bc applied to this set
      if( bcPresCompStatusMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Pressure boundary condition not prescribed on set {}", setName ) );
      }
    } );

    // 3. Check composition BC (global component fraction)
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    fields::flow::globalCompFraction::key(),
                                    [&] ( FieldSpecificationBase const & fs,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const &,
                                          FaceManager &,
                                          string const & )
    {
      // 3.1 Check pressure, temperature, and record composition bc application
      integer const comp = fs.getComponent();

      if( bcPresCompStatusMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Pressure boundary condition not prescribed on set {}", setName ) );
      }
      if( bcTempStatusMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Temperature boundary condition not prescribed on set {}. \n"
                                "Note that for face boundary conditions, you must provide a temperature", setName ) );
      }
      if( comp < 0 || comp >= m_numComponents )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Invalid component index [{}] in composition boundary condition {}", comp, fs.getName() ) );
        return; // can't check next part with invalid component id
      }

      ComponentMask< MAX_NC > & compMask = bcPresCompStatusMap[setName];
      if( compMask[comp] )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Conflicting composition[{}] boundary conditions on set {}", comp, setName ) );
      }
      compMask.set( comp );
    } );

    // 3.2 Check consistency between composition BC applied to sets
    for( auto const & setEntry : bcPresCompStatusMap )
    {
      ComponentMask< MAX_NC > const & compMask = setEntry.second;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        if( !compMask[ic] )
        {
          bcConsistent = false;
          GEOS_WARNING( GEOS_FMT( "Boundary condition not applied to composition[{}] on set {}", ic, setEntry.first ) );
        }
      }
    }
  } );

  return bcConsistent;
}

namespace
{
char const faceBcLogMessage[] =
  "CompositionalMultiphaseFVM {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target faces (including ghost faces) is {}."
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";
}

void CompositionalMultiphaseFVM::applyFaceDirichletBC( real64 const time_n,
                                                       real64 const dt,
                                                       DofManager const & dofManager,
                                                       DomainPartition & domain,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Only validate BC at the beginning of Newton loop
  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    bool const bcConsistent = validateFaceDirichletBC( domain, time_n + dt );
    GEOS_ERROR_IF( !bcConsistent, GEOS_FMT( "{}: inconsistent boundary conditions", getDataContext() ) );
  }

  using namespace isothermalCompositionalMultiphaseFVMKernels;

  // for now, we neglect capillary pressure in the kernel
  BitFlags< KernelFlags > kernelFlags;
  if( m_useTotalMassEquation )
    kernelFlags.set( KernelFlags::TotalMassEquation );

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // Take BCs defined for "pressure" field and apply values to "facePressure"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    fields::flow::pressure::key(), fields::flow::facePressure::key() );
    // Take BCs defined for "globalCompFraction" field and apply values to "faceGlobalCompFraction"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    fields::flow::globalCompFraction::key(), fields::flow::faceGlobalCompFraction::key() );
    // Take BCs defined for "temperature" field and apply values to "faceTemperature"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    fields::flow::temperature::key(), fields::flow::faceTemperature::key() );

    // Then launch the face Dirichlet kernel
    fsManager.apply< FaceManager >( time_n + dt,
                                    mesh,
                                    fields::flow::pressure::key(), // we have required that pressure is always present
                                    [&] ( FieldSpecificationBase const &,
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

      // TODO: same issue as in the single-phase case
      //       currently we just use model from the first cell in this stencil
      //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
      //       Can we just use cell properties for an approximate flux computation?
      //       Then we can forget about capturing the fluid model.
      localIndex const er = stencil.getElementRegionIndices()( 0, 0 );
      localIndex const esr = stencil.getElementSubRegionIndices()( 0, 0 );
      ElementSubRegionBase & subRegion = elemManager.getRegion( er ).getSubRegion( esr );
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & multiFluidBase = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      BoundaryStencilWrapper const stencilWrapper = stencil.createKernelWrapper();

      string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      if( m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
      {
        isothermalCompositionalMultiphaseFVMKernels::
          DirichletFluxComputeZFormulationKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     kernelFlags,
                                                     elemDofKey,
                                                     getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     multiFluidBase,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );

      }
      else
      {
        if( m_isThermal )
        {
          //todo (jafranc) extend upwindScheme name if satisfied in isothermalCase
          thermalCompositionalMultiphaseFVMKernels::
            DirichletFluxComputeKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       dofManager.rankOffset(),
                                                       kernelFlags,
                                                       elemDofKey,
                                                       getName(),
                                                       faceManager,
                                                       elemManager,
                                                       stencilWrapper,
                                                       multiFluidBase,
                                                       dt,
                                                       localMatrix,
                                                       localRhs );
        }
        else
        {
          isothermalCompositionalMultiphaseFVMKernels::
            DirichletFluxComputeKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       dofManager.rankOffset(),
                                                       kernelFlags,
                                                       elemDofKey,
                                                       getName(),
                                                       faceManager,
                                                       elemManager,
                                                       stencilWrapper,
                                                       multiFluidBase,
                                                       dt,
                                                       localMatrix,
                                                       localRhs );
        }
      }

    } );
  } );
}

void CompositionalMultiphaseFVM::applyAquiferBC( real64 const time,
                                                 real64 const dt,
                                                 DofManager const & dofManager,
                                                 DomainPartition & domain,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    isothermalCompositionalMultiphaseFVMKernels::
      AquiferBCKernel::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      AquiferBCKernel::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getName() );

    fsManager.apply< FaceManager,
                     AquiferBoundaryCondition >( time + dt,
                                                 mesh,
                                                 AquiferBoundaryCondition::catalogName(),
                                                 [&] ( AquiferBoundaryCondition const & bc,
                                                       string const & setName,
                                                       SortedArrayView< localIndex const > const &,
                                                       FaceManager & faceManager,
                                                       string const & )
    {
      BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
      if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
        GEOS_LOG_LEVEL_RANK_0_ON_GROUP( logInfo::BoundaryCondition,
                                        GEOS_FMT( faceBcLogMessage,
                                                  getName(), time+dt, bc.getCatalogName(), bc.getName(),
                                                  setName, faceManager.getName(), bc.getScale(), numTargetFaces ),
                                        bc );
      }

      if( stencil.size() == 0 )
      {
        return;
      }

      AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();
      bool const allowAllPhasesIntoAquifer = bc.allowAllPhasesIntoAquifer();
      localIndex const waterPhaseIndex = bc.getWaterPhaseIndex();
      real64 const & aquiferWaterPhaseDens = bc.getWaterPhaseDensity();
      arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac = bc.getWaterPhaseComponentFraction();

      // While this kernel is waiting for a factory class, pass all the accessors here
      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector1
      < isothermalCompositionalMultiphaseFVMKernels::AquiferBCKernel >( m_numComponents,
                                                                        m_numPhases,
                                                                        waterPhaseIndex,
                                                                        allowAllPhasesIntoAquifer,
                                                                        m_useTotalMassEquation,
                                                                        stencil,
                                                                        dofManager.rankOffset(),
                                                                        elemDofNumber.toNestedViewConst(),
                                                                        aquiferBCWrapper,
                                                                        aquiferWaterPhaseDens,
                                                                        aquiferWaterPhaseCompFrac,
                                                                        compFlowAccessors.get( fields::ghostRank{} ),
                                                                        compFlowAccessors.get( fields::flow::pressure{} ),
                                                                        compFlowAccessors.get( fields::flow::pressure_n{} ),
                                                                        compFlowAccessors.get( fields::flow::gravityCoefficient{} ),
                                                                        compFlowAccessors.get( fields::flow::phaseVolumeFraction{} ),
                                                                        compFlowAccessors.get( fields::flow::dPhaseVolumeFraction{} ),
                                                                        compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity{} ),
                                                                        multiFluidAccessors.get( fields::multifluid::phaseDensity{} ),
                                                                        multiFluidAccessors.get( fields::multifluid::dPhaseDensity{} ),
                                                                        multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ),
                                                                        multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction{} ),
                                                                        time,
                                                                        dt,
                                                                        localMatrix.toViewConstSizes(),
                                                                        localRhs.toView() );
    } );
  } );

}

void CompositionalMultiphaseFVM::assembleHydrofracFluxTerms( real64 const GEOS_UNUSED_PARAM ( time_n ),
                                                             real64 const dt,
                                                             DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs,
                                                             CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  using namespace isothermalCompositionalMultiphaseFVMKernels;

  BitFlags< KernelFlags > kernelFlags;
  if( m_hasCapPressure )
    kernelFlags.set( KernelFlags::CapPressure );
  if( m_useTotalMassEquation )
    kernelFlags.set( KernelFlags::TotalMassEquation );
  if( m_gravityDensityScheme == GravityDensityScheme::PhasePresence )
    kernelFlags.set( KernelFlags::CheckPhasePresenceInGravity );

  if( m_isThermal )
  {
    // should not end up here but just in case
    GEOS_ERROR( "Thermal not yet supported in CompositionalMultiphaseFVM::assembleHydrofracFluxTerms" );
  }
  if( fluxApprox.upwindingParams().upwindingScheme != UpwindingScheme::PPU )
  {
    // a bit tricky to check in advance
    GEOS_ERROR( "Only PPU upwinding is supported in CompositionalMultiphaseFVM::assembleHydrofracFluxTerms" );
  }

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & )
  {
    fluxApprox.forStencils< CellElementStencilTPFA, FaceElementToCellStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      isothermalCompositionalMultiphaseFVMKernels::
        FluxComputeKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   m_numPhases,
                                                   dofManager.rankOffset(),
                                                   elemDofKey,
                                                   kernelFlags,
                                                   getName(),
                                                   mesh.getElemManager(),
                                                   stencilWrapper,
                                                   dt,
                                                   localMatrix.toViewConstSizes(),
                                                   localRhs.toView() );
    } );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      multiphasePoromechanicsConformingFracturesKernels::
        FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                                             m_numPhases,
                                                                             dofManager.rankOffset(),
                                                                             elemDofKey,
                                                                             kernelFlags,
                                                                             getName(),
                                                                             mesh.getElemManager(),
                                                                             stencilWrapper,
                                                                             dt,
                                                                             localMatrix.toViewConstSizes(),
                                                                             localRhs.toView(),
                                                                             dR_dAper );
    } );
  } );
}

real64 CompositionalMultiphaseFVM::setNextDt( real64 const & currentTime,
                                              real64 const & currentDt,
                                              DomainPartition & domain )
{
  if( m_targetFlowCFL < 0 )
    return CompositionalMultiphaseBase::setNextDt( currentTime, currentDt, domain );
  else
    return setNextDtBasedOnCFL( currentDt, domain );
}

real64 CompositionalMultiphaseFVM::setNextDtBasedOnCFL( const geos::real64 & currentDt, geos::DomainPartition & domain )
{

  real64 maxPhaseCFL, maxCompCFL;

  computeCFLNumbers( domain, currentDt, maxPhaseCFL, maxCompCFL );

  GEOS_LOG_LEVEL_RANK_0( logInfo::TimeStep,
                         GEOS_FMT( "{}: max phase CFL number = {}", getName(), maxPhaseCFL ) );
  GEOS_LOG_LEVEL_RANK_0( logInfo::TimeStep,
                         GEOS_FMT( "{}: max component CFL number = {} ", getName(), maxCompCFL ) );

  return std::min( m_targetFlowCFL * currentDt / maxCompCFL,
                   m_targetFlowCFL * currentDt / maxPhaseCFL );
}

void CompositionalMultiphaseFVM::computeCFLNumbers( geos::DomainPartition & domain, const geos::real64 & dt,
                                                    geos::real64 & maxPhaseCFL, geos::real64 & maxCompCFL )
{
  GEOS_MARK_FUNCTION;

  integer const numPhases = numFluidPhases();
  integer const numComps = numFluidComponents();

  // Step 1: reset the arrays involved in the computation of CFL numbers
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView2d< real64, compflow::USD_PHASE > const & phaseOutflux =
        subRegion.getField< fields::flow::phaseOutflux >();
      arrayView2d< real64, compflow::USD_COMP > const & compOutflux =
        subRegion.getField< fields::flow::componentOutflux >();
      phaseOutflux.zero();
      compOutflux.zero();
    } );

    // Step 2: compute the total volumetric outflux for each reservoir cell by looping over faces
    NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( getDiscretizationName() );

    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::PermeabilityAccessors permeabilityAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::RelPermAccessors relPermAccessors( mesh.getElemManager(), getName() );

    // TODO: find a way to compile with this modifiable accessors in CompFlowAccessors, and remove them from here
    ElementRegionManager::ElementViewAccessor< arrayView2d< real64, compflow::USD_PHASE > > const phaseOutfluxAccessor =
      mesh.getElemManager().constructViewAccessor< array2d< real64, compflow::LAYOUT_PHASE >,
                                                   arrayView2d< real64, compflow::USD_PHASE > >( fields::flow::phaseOutflux::key() );

    ElementRegionManager::ElementViewAccessor< arrayView2d< real64, compflow::USD_COMP > > const compOutfluxAccessor =
      mesh.getElemManager().constructViewAccessor< array2d< real64, compflow::LAYOUT_COMP >,
                                                   arrayView2d< real64, compflow::USD_COMP > >( fields::flow::componentOutflux::key() );


    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      // While this kernel is waiting for a factory class, pass all the accessors here
      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector1
      < isothermalCompositionalMultiphaseFVMKernels::CFLFluxKernel >( numComps,
                                                                      numPhases,
                                                                      m_gravityDensityScheme == GravityDensityScheme::PhasePresence,
                                                                      dt,
                                                                      stencilWrapper,
                                                                      compFlowAccessors.get( fields::flow::pressure{} ),
                                                                      compFlowAccessors.get( fields::flow::gravityCoefficient{} ),
                                                                      compFlowAccessors.get( fields::flow::phaseVolumeFraction{} ),
                                                                      permeabilityAccessors.get( fields::permeability::permeability{} ),
                                                                      permeabilityAccessors.get( fields::permeability::dPerm_dPressure{} ),
                                                                      relPermAccessors.get( fields::relperm::phaseRelPerm{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseViscosity{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseDensity{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseMassDensity{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ),
                                                                      phaseOutfluxAccessor.toNestedView(),
                                                                      compOutfluxAccessor.toNestedView() );
    } );
  } );

  // Step 3: finalize the (cell-based) computation of the CFL numbers
  real64 localMaxPhaseCFLNumber = 0.0;
  real64 localMaxCompCFLNumber = 0.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux =
        subRegion.getField< fields::flow::phaseOutflux >();
      arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux =
        subRegion.getField< fields::flow::componentOutflux >();

      arrayView1d< real64 > const & phaseCFLNumber = subRegion.getField< fields::flow::phaseCFLNumber >();
      arrayView1d< real64 > const & compCFLNumber = subRegion.getField< fields::flow::componentCFLNumber >();

      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

      arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
        subRegion.getField< fields::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
        subRegion.getField< fields::flow::globalCompFraction >();
      arrayView2d< real64, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::flow::phaseVolumeFraction >();

      Group const & constitutiveModels = subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

      string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = constitutiveModels.getGroup< MultiFluidBase >( fluidName );
      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc = fluid.phaseViscosity();

      string const & relpermName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relperm = constitutiveModels.getGroup< RelativePermeabilityBase >( relpermName );
      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm = relperm.phaseRelPerm();
      arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac = relperm.dPhaseRelPerm_dPhaseVolFraction();

      string const & solidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = constitutiveModels.getGroup< CoupledSolidBase >( solidName );
      arrayView2d< real64 const > const & porosity = solid.getPorosity();

      real64 subRegionMaxPhaseCFLNumber = 0.0;
      real64 subRegionMaxCompCFLNumber = 0.0;

      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector2
      < isothermalCompositionalMultiphaseFVMKernels::CFLKernel >( numComps, numPhases,
                                                                  subRegion.size(),
                                                                  volume,
                                                                  porosity,
                                                                  compDens,
                                                                  compFrac,
                                                                  phaseVolFrac,
                                                                  phaseRelPerm,
                                                                  dPhaseRelPerm_dPhaseVolFrac,
                                                                  phaseVisc,
                                                                  phaseOutflux,
                                                                  compOutflux,
                                                                  phaseCFLNumber,
                                                                  compCFLNumber,
                                                                  subRegionMaxPhaseCFLNumber,
                                                                  subRegionMaxCompCFLNumber );

      localMaxPhaseCFLNumber = LvArray::math::max( localMaxPhaseCFLNumber, subRegionMaxPhaseCFLNumber );
      localMaxCompCFLNumber = LvArray::math::max( localMaxCompCFLNumber, subRegionMaxCompCFLNumber );

    } );
  } );

  maxPhaseCFL = MpiWrapper::max( localMaxPhaseCFLNumber );
  maxCompCFL = MpiWrapper::max( localMaxCompCFLNumber );

}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, CompositionalMultiphaseFVM, string const &, Group * const )
//END_SPHINX_INCLUDE_01
}// namespace geos
