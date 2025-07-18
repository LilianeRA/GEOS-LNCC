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
 * @file SinglePhaseBase.cpp
 */

#include "SinglePhaseBase.hpp"


#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidSelector.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/solid/SolidInternalEnergy.hpp"
#include "constitutive/thermalConductivity/SinglePhaseThermalConductivitySelector.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/LogLevelsInfo.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "physicsSolvers/fluidFlow/SourceFluxStatistics.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "functions/TableFunction.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/KernelLaunchSelectors.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/MobilityKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/SolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/SolutionScalingKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/StatisticsKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/HydrostaticPressureKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluidUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/SolidInternalEnergyUpdateKernel.hpp"


namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace singlePhaseBaseKernels;

SinglePhaseBase::SinglePhaseBase( const string & name,
                                  Group * const parent ):
  FlowSolverBase( name, parent )
{
  this->registerWrapper( viewKeyStruct::inputTemperatureString(), &m_inputTemperature ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature" );

  this->getWrapper< integer >( string( viewKeyStruct::isThermalString() ) ).
    appendDescription( GEOS_FMT( "\nSourceFluxes application if {} is enabled :\n"
                                 "- negative value (injection): the mass balance equation is modified to considered the additional source term,\n"
                                 "- positive value (production): both the mass balance and the energy balance equations are modified to considered the additional source term.\n"
                                 "For the energy balance equation, the mass flux is multiplied by the enthalpy in the cell from which the fluid is being produced.",
                                 viewKeyStruct::isThermalString() ) );

  addLogLevel< logInfo::Solution >();
  addLogLevel< logInfo::SourceFluxFailure >();
}


void SinglePhaseBase::registerDataOnMesh( Group & meshBodies )
{
  using namespace fields::flow;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  m_numDofPerCell = m_isThermal ? 2 : 1;

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< mobility >( getName() );
      subRegion.registerField< dMobility >( getName()).reference().resizeDimension< 1 >( m_numDofPerCell );

      subRegion.registerField< mass >( getName() );
      subRegion.registerField< mass_n >( getName() );
      subRegion.registerField< dMass >( getName() ).reference().resizeDimension< 1 >( m_numDofPerCell );

      if( m_isThermal )
      {
        subRegion.registerField< dEnergy >( getName() ).reference().resizeDimension< 1 >( m_numDofPerCell );
      }
    } );

    elemManager.forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      SurfaceElementSubRegion & subRegion )
    {
      subRegion.registerField< massCreated >( getName() );
    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      faceManager.registerField< facePressure >( getName() );

      if( m_isThermal )
      {
        faceManager.registerField< faceTemperature >( getName() );
      }
    }
  } );
}

void SinglePhaseBase::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  setConstitutiveName< SingleFluidBase >( subRegion, viewKeyStruct::fluidNamesString(), "singlephase fluid" );

  if( m_isThermal )
  {
    setConstitutiveName< SinglePhaseThermalConductivityBase >( subRegion, viewKeyStruct::thermalConductivityNamesString(), "singlephase thermal conductivity" );
  }
}

void SinglePhaseBase::initializeAquiferBC() const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition & bc )
  {
    // set the gravity vector (needed later for the potential diff calculations)
    bc.setGravityVector( gravityVector() );
  } );
}

void SinglePhaseBase::validateConstitutiveModels( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      fluidName = getConstitutiveName< SingleFluidBase >( subRegion );
      GEOS_THROW_IF( fluidName.empty(),
                     GEOS_FMT( "SingleFluidBase {}: Fluid model not found on subregion {}",
                               getDataContext(), subRegion.getName() ),
                     InputError );

      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

      constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
      {
        string const fluidModelName = castedFluid.getCatalogName();
        GEOS_THROW_IF( m_isThermal && (fluidModelName != "ThermalCompressibleSinglePhaseFluid"),
                       GEOS_FMT( "SingleFluidBase {}: the thermal option is enabled in the solver, but the fluid model {} is not for thermal fluid",
                                 getDataContext(), fluid.getDataContext() ),
                       InputError );
        GEOS_THROW_IF( !m_isThermal && (fluidModelName == "ThermalCompressibleSinglePhaseFluid"),
                       GEOS_FMT( "SingleFluidBase {}: the fluid model is for thermal fluid {}, but the solver option is incompatible with the fluid model",
                                 getDataContext(), fluid.getDataContext() ),
                       InputError );
      } );
    } );
  } );
}

void SinglePhaseBase::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // 1. Validate various models against each other (must have same phases and components)
  validateConstitutiveModels( domain );

  // 2. Set the value of temperature
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const temp = subRegion.getField< fields::flow::temperature >();
      temp.setValues< parallelHostPolicy >( m_inputTemperature );
    } );
  } );

  // 3. Initialize the aquifer boundary condition
  initializeAquiferBC();
}

void SinglePhaseBase::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const temp = dataGroup.getField< fields::flow::temperature >();

  SingleFluidBase & fluid =
    getConstitutiveModel< SingleFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    singlePhaseBaseKernels::FluidUpdateKernel::launch( fluidWrapper, pres, temp );
  } );
}

void SinglePhaseBase::updateMass( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  geos::internal::kernelLaunchSelectorThermalSwitch ( m_isThermal, [&] ( auto ISTHERMAL )
  {
    integer constexpr IS_THERMAL = ISTHERMAL();
    updateMass< IS_THERMAL >( subRegion );
  } );
}

template< integer IS_THERMAL >
void SinglePhaseBase::updateMass( ElementSubRegionBase & subRegion ) const
{
  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

  arrayView1d< real64 > const mass = subRegion.getField< fields::flow::mass >();
  arrayView1d< real64 > const mass_n = subRegion.getField< fields::flow::mass_n >();
  arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const dMass = subRegion.getField< fields::flow::dMass >();

  //START_SPHINX_INCLUDE_COUPLEDSOLID
  CoupledSolidBase const & porousSolid =
    getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString()));
  //END_SPHINX_INCLUDE_COUPLEDSOLID
  arrayView2d< real64 const > const porosity = porousSolid.getPorosity();
  arrayView2d< real64 const > const dPorosity_dP = porousSolid.getDporosity_dPressure();
  arrayView2d< real64 const > const porosity_n = porousSolid.getPorosity_n();

  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView1d< real64 const > const deltaVolume = subRegion.getField< fields::flow::deltaVolume >();

  SingleFluidBase & fluid =
    getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString()));
  arrayView2d< real64 const, singlefluid::USD_FLUID > const density = fluid.density();
  arrayView2d< real64 const, singlefluid::USD_FLUID > const density_n = fluid.density_n();
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dDensity = fluid.dDensity();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 const vol = volume[ei] + deltaVolume[ei];
    mass[ei] = porosity[ei][0] * density[ei][0] * vol;
    dMass[ei][DerivOffset::dP] = ( dPorosity_dP[ei][0] * density[ei][0] + porosity[ei][0] * dDensity[ei][0][DerivOffset::dP] ) * vol;
    if( isZero( mass_n[ei] ) )   // this is a hack for hydrofrac cases
    {
      mass_n[ei] = porosity_n[ei][0] * volume[ei] * density_n[ei][0];   // initialize newly created element mass
    }
  } );

  if constexpr (IS_THERMAL)
  {
    arrayView2d< real64 const > const dPorosity_dT = porousSolid.getDporosity_dTemperature();
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      real64 const vol = volume[ei] + deltaVolume[ei];
      dMass[ei][DerivOffset::dT] = ( dPorosity_dT[ei][0] * density[ei][0] + porosity[ei][0] * dDensity[ei][0][DerivOffset::dT] ) * vol;
    } );
  }
}

void SinglePhaseBase::updateEnergy( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;

  arrayView1d< real64 > const energy = subRegion.getField< fields::flow::energy >();
  arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const dEnergy = subRegion.getField< fields::flow::dEnergy >();

  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView1d< real64 const > const deltaVolume = subRegion.getField< fields::flow::deltaVolume >();

  CoupledSolidBase const & porousSolid =
    getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
  arrayView2d< real64 const > const porosity = porousSolid.getPorosity();
  arrayView2d< real64 const > const dPorosity_dP = porousSolid.getDporosity_dPressure();
  arrayView2d< real64 const > const dPorosity_dT = porousSolid.getDporosity_dTemperature();
  arrayView2d< real64 const > const rockInternalEnergy = porousSolid.getInternalEnergy();
  arrayView2d< real64 const > const dRockInternalEnergy_dT = porousSolid.getDinternalEnergy_dTemperature();

  SingleFluidBase & fluid =
    getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const density = fluid.density();
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const fluidInternalEnergy = fluid.internalEnergy();
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dDensity = fluid.dDensity();
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dFluidInternalEnergy = fluid.dInternalEnergy();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 const vol = volume[ei] + deltaVolume[ei];
    energy[ei] = vol *
                 ( porosity[ei][0] * density[ei][0] * fluidInternalEnergy[ei][0] +
                   ( 1.0 - porosity[ei][0] ) * rockInternalEnergy[ei][0] );
    dEnergy[ei][DerivOffset::dP] = vol *
                                   ( dPorosity_dP[ei][0] * density[ei][0] * fluidInternalEnergy[ei][0] +
                                     porosity[ei][0] * dDensity[ei][0][DerivOffset::dP] * fluidInternalEnergy[ei][0] +
                                     porosity[ei][0] * density[ei][0] * dFluidInternalEnergy[ei][0][DerivOffset::dP] -
                                     dPorosity_dP[ei][0] * rockInternalEnergy[ei][0] );
    dEnergy[ei][DerivOffset::dT] = vol *
                                   ( dPorosity_dT[ei][0] * density[ei][0] * fluidInternalEnergy[ei][0] +
                                     porosity[ei][0] * dDensity[ei][0][DerivOffset::dT] * fluidInternalEnergy[ei][0] +
                                     porosity[ei][0] * density[ei][0] * dFluidInternalEnergy[ei][0][DerivOffset::dT] -
                                     dPorosity_dT[ei][0] * rockInternalEnergy[ei][0] +
                                     ( 1.0 - porosity[ei][0] ) * dRockInternalEnergy_dT[ei][0] );
  } );
}

void SinglePhaseBase::updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const
{
  arrayView1d< real64 const > const temperature = dataGroup.getField< fields::flow::temperature >();

  string const & solidInternalEnergyName = dataGroup.getReference< string >( viewKeyStruct::solidInternalEnergyNamesString() );
  SolidInternalEnergy & solidInternalEnergy = getConstitutiveModel< SolidInternalEnergy >( dataGroup, solidInternalEnergyName );

  SolidInternalEnergy::KernelWrapper solidInternalEnergyWrapper = solidInternalEnergy.createKernelUpdates();

  thermalSinglePhaseBaseKernels::SolidInternalEnergyUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(), solidInternalEnergyWrapper, temperature );
}

void SinglePhaseBase::updateThermalConductivity( ElementSubRegionBase & subRegion ) const
{
  string const & thermalConductivityName = subRegion.template getReference< string >( viewKeyStruct::thermalConductivityNamesString() );
  SinglePhaseThermalConductivityBase const & conductivityMaterial =
    getConstitutiveModel< SinglePhaseThermalConductivityBase >( subRegion, thermalConductivityName );

  arrayView1d< real64 const > const & temperature = subRegion.template getField< fields::flow::temperature >();
  conductivityMaterial.updateFromTemperature( temperature );
}

real64 SinglePhaseBase::updateFluidState( ElementSubRegionBase & subRegion ) const
{
  updateFluidModel( subRegion );
  updateMass( subRegion );
  updateMobility( subRegion );
  return 0.0;
}

void SinglePhaseBase::updateMobility( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  // output
  arrayView1d< real64 > const mob = dataGroup.getField< fields::flow::mobility >();
  arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const dMobility = dataGroup.getField< fields::flow::dMobility >();

  // input
  SingleFluidBase & fluid =
    getConstitutiveModel< SingleFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  geos::internal::kernelLaunchSelectorThermalSwitch( m_isThermal, [&] ( auto ISTHERMAL ) {
    integer constexpr NUMDOF = ISTHERMAL() + 1;
    singlePhaseBaseKernels::MobilityKernel::compute_value_and_derivatives< parallelDevicePolicy<>, NUMDOF >( dataGroup.size(),
                                                                                                             fluid.density(),
                                                                                                             fluid.dDensity(),
                                                                                                             fluid.viscosity(),
                                                                                                             fluid.dViscosity(),
                                                                                                             mob,
                                                                                                             dMobility );
  } );

}

void SinglePhaseBase::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  allowNegativePressure();

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  FlowSolverBase::initializeState( domain );
}

void SinglePhaseBase::computeHydrostaticEquilibrium( DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  // Step 1: count individual equilibriums (there may be multiple ones)

  std::map< string, localIndex > equilNameToEquilId;
  localIndex equilCounter = 0;

  fsManager.forSubGroups< EquilibriumInitialCondition >( [&] ( EquilibriumInitialCondition const & bc )
  {
    // collect all the equil name to idx
    equilNameToEquilId[bc.getName()] = equilCounter;
    equilCounter++;

    // check that the gravity vector is aligned with the z-axis
    GEOS_THROW_IF( !isZero( gravVector[0] ) || !isZero( gravVector[1] ),
                   getCatalogName() << " " << getDataContext() <<
                   ": the gravity vector specified in this simulation (" << gravVector[0] << " " << gravVector[1] << " " << gravVector[2] <<
                   ") is not aligned with the z-axis. \n"
                   "This is incompatible with the " << bc.getCatalogName() << " " << bc.getDataContext() <<
                   "used in this simulation. To proceed, you can either: \n" <<
                   "   - Use a gravityVector aligned with the z-axis, such as (0.0,0.0,-9.81)\n" <<
                   "   - Remove the hydrostatic equilibrium initial condition from the XML file",
                   InputError );
  } );

  if( equilCounter == 0 )
  {
    return;
  }

  // Step 2: find the min elevation and the max elevation in the targetSets
  array1d< real64 > globalMaxElevation( equilNameToEquilId.size() );
  array1d< real64 > globalMinElevation( equilNameToEquilId.size() );
  findMinMaxElevationInEquilibriumTarget( domain,
                                          equilNameToEquilId,
                                          globalMaxElevation,
                                          globalMinElevation );

  // Step 3: for each equil, compute a fine table with hydrostatic pressure vs elevation if the region is a target region
  // first compute the region filter
  std::set< string > regionFilter;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel &,
                                                                string_array const & regionNames )
  {
    for( string const & regionName : regionNames )
    {
      regionFilter.insert( regionName );
    }
  } );

  // then start the actual table construction
  fsManager.apply< ElementSubRegionBase,
                   EquilibriumInitialCondition >( 0.0,
                                                  domain.getMeshBody( 0 ).getBaseDiscretization(),
                                                  EquilibriumInitialCondition::catalogName(),
                                                  [&] ( EquilibriumInitialCondition const & fs,
                                                        string const &,
                                                        SortedArrayView< localIndex const > const & targetSet,
                                                        ElementSubRegionBase & subRegion,
                                                        string const & )
  {
    // Step 3.1: retrieve the data necessary to construct the pressure table in this subregion

    integer const maxNumEquilIterations = fs.getMaxNumEquilibrationIterations();
    real64 const equilTolerance = fs.getEquilibrationTolerance();
    real64 const datumElevation = fs.getDatumElevation();
    real64 const datumPressure = fs.getDatumPressure();

    localIndex const equilIndex = equilNameToEquilId.at( fs.getName() );
    real64 const minElevation = LvArray::math::min( globalMinElevation[equilIndex], datumElevation );
    real64 const maxElevation = LvArray::math::max( globalMaxElevation[equilIndex], datumElevation );
    real64 const elevationIncrement = LvArray::math::min( fs.getElevationIncrement(), maxElevation - minElevation );
    localIndex const numPointsInTable = ( elevationIncrement > 0 ) ? std::ceil( (maxElevation - minElevation) / elevationIncrement ) + 1 : 1;

    real64 const eps = 0.1 * (maxElevation - minElevation); // we add a small buffer to only log in the pathological cases
    GEOS_LOG_RANK_0_IF( ( (datumElevation > globalMaxElevation[equilIndex]+eps)  || (datumElevation < globalMinElevation[equilIndex]-eps) ),
                        getCatalogName() << " " << getDataContext() <<
                        ": By looking at the elevation of the cell centers in this model, GEOS found that " <<
                        "the min elevation is " << globalMinElevation[equilIndex] << " and the max elevation is " <<
                        globalMaxElevation[equilIndex] << "\nBut, a datum elevation of " << datumElevation <<
                        " was specified in the input file to equilibrate the model.\n " <<
                        "The simulation is going to proceed with this out-of-bound datum elevation," <<
                        " but the initial condition may be inaccurate." );

    array1d< array1d< real64 > > elevationValues;
    array1d< real64 > pressureValues;
    elevationValues.resize( 1 );
    elevationValues[0].resize( numPointsInTable );
    pressureValues.resize( numPointsInTable );

    // Step 3.2: retrieve the fluid model to compute densities
    // we end up with the same issue as in applyDirichletBC: there is not a clean way to retrieve the fluid info

    // filter out region not in target
    Group const & region = subRegion.getParent().getParent();
    auto it = regionFilter.find( region.getName() );
    if( it == regionFilter.end() )
    {
      return; // the region is not in target, there is nothing to do
    }

    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());

    // filter out the proppant fluid constitutive models
    ConstitutiveBase & fluid = getConstitutiveModel( subRegion, fluidName );
    if( !dynamicCast< SingleFluidBase * >( &fluid ) )
    {
      return;
    }
    SingleFluidBase & singleFluid = dynamicCast< SingleFluidBase & >( fluid );

    // Step 3.3: compute the hydrostatic pressure values

    constitutiveUpdatePassThru( singleFluid, [&] ( auto & castedFluid )
    {
      using FluidType = TYPEOFREF( castedFluid );
      typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      // note: inside this kernel, serialPolicy is used, and elevation/pressure values don't go to the GPU
      bool const equilHasConverged =
        HydrostaticPressureKernel::launch( numPointsInTable,
                                           maxNumEquilIterations,
                                           equilTolerance,
                                           gravVector,
                                           minElevation,
                                           elevationIncrement,
                                           datumElevation,
                                           datumPressure,
                                           fluidWrapper,
                                           elevationValues.toNestedView(),
                                           pressureValues.toView() );

      GEOS_THROW_IF( !equilHasConverged,
                     getCatalogName() << " " << getDataContext() <<
                     ": hydrostatic pressure initialization failed to converge in region " << region.getName() << "!",
                     std::runtime_error );
    } );

    // Step 3.4: create hydrostatic pressure table

    FunctionManager & functionManager = FunctionManager::getInstance();

    string const tableName = fs.getName() + "_" + subRegion.getName() + "_table";
    TableFunction * const presTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
    presTable->setTableCoordinates( elevationValues, { units::Distance } );
    presTable->setTableValues( pressureValues, units::Pressure );
    presTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    TableFunction::KernelWrapper presTableWrapper = presTable->createKernelWrapper();

    // Step 4: assign pressure as a function of elevation
    // TODO: this last step should probably be delayed to wait for the creation of FaceElements
    arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();
    arrayView1d< real64 > const pres = subRegion.getField< fields::flow::pressure >();

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minPressure( LvArray::NumericLimits< real64 >::max );

    forAll< parallelDevicePolicy< > >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const k = targetSet[i];
      real64 const elevation = elemCenter[k][2];
      pres[k] = presTableWrapper.compute( &elevation );
      minPressure.min( pres[k] );
    } );

    // For single phase flow, just issue a warning, because the simulation can proceed with a negative pressure
    GEOS_WARNING_IF( minPressure.get() <= 0.0,
                     GEOS_FMT( "A negative pressure of {} Pa was found during hydrostatic initialization in region/subRegion {}/{}",
                               minPressure.get(), region.getName(), subRegion.getName() ) );
  } );
}

void SinglePhaseBase::initializeFluidState( MeshLevel & mesh, string_array const & regionNames )
{
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                 auto & subRegion )
  {
    SingleFluidBase const & fluid =
      getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString()));
    updateFluidState( subRegion );

    // 2. save the initial density (for use in the single-phase poromechanics solver to compute the deltaBodyForce)
    fluid.initializeState();
  } );
}

void SinglePhaseBase::initializeThermalState( MeshLevel & mesh, string_array const & regionNames )
{
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                 auto & subRegion )
  {
    // initialized porosity
    CoupledSolidBase const & porousSolid =
      getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
    arrayView2d< real64 const > const porosity = porousSolid.getPorosity();

    string const & thermalConductivityName = subRegion.template getReference< string >( viewKeyStruct::thermalConductivityNamesString());
    SinglePhaseThermalConductivityBase const & conductivityMaterial =
      getConstitutiveModel< SinglePhaseThermalConductivityBase >( subRegion, thermalConductivityName );
    conductivityMaterial.initializeRockFluidState( porosity );
    // note that there is nothing to update here because thermal conductivity is explicit for now

    updateSolidInternalEnergyModel( subRegion );
    string const & solidInternalEnergyName = subRegion.template getReference< string >( viewKeyStruct::solidInternalEnergyNamesString());
    SolidInternalEnergy const & solidInternalEnergyMaterial =
      getConstitutiveModel< SolidInternalEnergy >( subRegion, solidInternalEnergyName );
    solidInternalEnergyMaterial.saveConvergedState();

    updateEnergy( subRegion );
  } );
}

void SinglePhaseBase::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                         real64 const & GEOS_UNUSED_PARAM( dt ),
                                         DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      saveConvergedState( subRegion );

      applyDeltaVolume( subRegion );

      // This should fix NaN density in newly created fracture elements
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );
      // for thermal simulations, update solid internal energy
      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
        updateThermalConductivity( subRegion );
        updateEnergy( subRegion );
      }

    } );

    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                             SurfaceElementSubRegion & subRegion )
    {
      arrayView1d< real64 const > const aper = subRegion.getField< fields::flow::hydraulicAperture >();
      arrayView1d< real64 > const aper0 = subRegion.getField< fields::flow::aperture0 >();
      aper0.setValues< parallelDevicePolicy<> >( aper );

      // Needed coz faceElems don't exist when initializing.
      CoupledSolidBase const & porousSolid =
        getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::solidNamesString() ) );
      porousSolid.saveConvergedState();

      saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );

      // This call is required by the proppant solver, but should not be here
      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
      fluid.saveConvergedState();

    } );

  } );
}

void SinglePhaseBase::implicitStepComplete( real64 const & time,
                                            real64 const & dt,
                                            DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // note: we have to save the aquifer state **before** updating the pressure,
  // otherwise the aquifer flux is saved with the wrong pressure time level
  saveAquiferConvergedState( time, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      // update deltaPressure
      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const initPres = subRegion.getField< fields::flow::initialPressure >();
      arrayView1d< real64 > const deltaPres = subRegion.getField< fields::flow::deltaPressure >();
      singlePhaseBaseKernels::StatisticsKernel::
        saveDeltaPressure( subRegion.size(), pres, initPres, deltaPres );

      applyDeltaVolume( subRegion );

      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );
      fluid.saveConvergedState();

      CoupledSolidBase const & porousSolid =
        getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
      if( m_keepVariablesConstantDuringInitStep )
      {
        porousSolid.ignoreConvergedState(); // newPorosity <- porosity_n
      }
      else
      {
        porousSolid.saveConvergedState(); // porosity_n <- porosity
      }

      if( m_isThermal )
      {
        arrayView2d< real64 const > const porosity = porousSolid.getPorosity();

        SinglePhaseThermalConductivityBase const & conductivityMaterial =
          getConstitutiveModel< SinglePhaseThermalConductivityBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::thermalConductivityNamesString() ) );

        conductivityMaterial.saveConvergedRockFluidState( porosity );
      }

    } );

    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                             SurfaceElementSubRegion & subRegion )
    {
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();
      arrayView1d< real64 > const creationMass = subRegion.getField< fields::flow::massCreated >();

      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );
      arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const density_n = fluid.density_n();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          if( volume[ei] * density_n[ei][0] > 1.1 * creationMass[ei] )
          {
            creationMass[ei] *= 0.75;
            if( creationMass[ei] < 1.0e-20 )
            {
              creationMass[ei] = 0.0;
            }
          }
        }
      } );
    } );
  } );
}


void SinglePhaseBase::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                      real64 const dt,
                                      DomainPartition & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  assembleAccumulationTerms( domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  if( m_isJumpStabilized )
  {
    assembleStabilizedFluxTerms( dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );
  }
  else
  {
    assembleFluxTerms( dt,
                       domain,
                       dofManager,
                       localMatrix,
                       localRhs );
  }
}

void SinglePhaseBase::assembleAccumulationTerms( DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      accumulationAssemblyLaunch( dofManager, subRegion, localMatrix, localRhs );
    } );
  } );
}

void SinglePhaseBase::applyBoundaryConditions( real64 time_n,
                                               real64 dt,
                                               DomainPartition & domain,
                                               DofManager const & dofManager,
                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                               arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  if( m_keepVariablesConstantDuringInitStep )
  {
    // this function is going to force the current flow state to be constant during the time step
    // this is used when the poromechanics solver is performing the stress initialization
    // TODO: in the future, a dedicated poromechanics kernel should eliminate the flow vars to construct a reduced system
    //       which will remove the need for this brittle passing aroung of flag
    keepVariablesConstantDuringInitStep( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
  }
  else
  {
    applySourceFluxBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
    applyDirichletBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
    applyAquiferBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }
}

namespace
{

void applyAndSpecifyFieldValue( real64 const & time_n,
                                real64 const & dt,
                                MeshLevel & mesh,
                                globalIndex const rankOffset,
                                string const dofKey,
                                bool const,
                                integer const idof,
                                string const fieldKey,
                                string const boundaryFieldKey,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                           mesh,
                                           fieldKey,
                                           [&]( FieldSpecificationBase const & fs,
                                                string const &,
                                                SortedArrayView< localIndex const > const & lset,
                                                ElementSubRegionBase & subRegion,
                                                string const & )
  {
    // Specify the bc value of the field
    fs.applyFieldValue< FieldSpecificationEqual,
                        parallelDevicePolicy<> >( lset,
                                                  time_n + dt,
                                                  subRegion,
                                                  boundaryFieldKey );

    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
    arrayView1d< globalIndex const > const dofNumber =
      subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< real64 const > const bcField =
      subRegion.getReference< array1d< real64 > >( boundaryFieldKey );
    arrayView1d< real64 const > const field =
      subRegion.getReference< array1d< real64 > >( fieldKey );

    forAll< parallelDevicePolicy<> >( lset.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      localIndex const ei = lset[a];
      if( ghostRank[ei] >= 0 )
      {
        return;
      }

      globalIndex const dofIndex = dofNumber[ei];
      localIndex const localRow = dofIndex - rankOffset;
      real64 rhsValue;

      // Apply field value to the matrix/rhs
      FieldSpecificationEqual::SpecifyFieldValue( dofIndex + idof,
                                                  rankOffset,
                                                  localMatrix,
                                                  rhsValue,
                                                  bcField[ei],
                                                  field[ei] );
      localRhs[localRow + idof] = rhsValue;
    } );
  } );
}

}

void SinglePhaseBase::applyDirichletBC( real64 const time_n,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  globalIndex const rankOffset = dofManager.rankOffset();
  bool const isFirstNonlinearIteration = ( m_nonlinearSolverParameters.m_numNewtonIterations == 0 );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    applyAndSpecifyFieldValue( time_n, dt, mesh, rankOffset, dofKey, isFirstNonlinearIteration,
                               0, fields::flow::pressure::key(), fields::flow::bcPressure::key(),
                               localMatrix, localRhs );
    if( m_isThermal )
    {
      applyAndSpecifyFieldValue( time_n, dt, mesh, rankOffset, dofKey, isFirstNonlinearIteration,
                                 1, fields::flow::temperature::key(), fields::flow::bcTemperature::key(),
                                 localMatrix, localRhs );
    }
  } );
}

void SinglePhaseBase::applySourceFluxBC( real64 const time_n,
                                         real64 const dt,
                                         DomainPartition & domain,
                                         DofManager const & dofManager,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Step 1: count individual source flux boundary conditions

  std::map< string, localIndex > bcNameToBcId;
  localIndex bcCounter = 0;

  fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&] ( SourceFluxBoundaryCondition const & bc )
  {
    // collect all the bc names to idx
    bcNameToBcId[bc.getName()] = bcCounter;
    bcCounter++;
  } );

  if( bcCounter == 0 )
  {
    return;
  }

  // Step 2: count the set size for each source flux (each source flux may have multiple target sets)

  array1d< globalIndex > bcAllSetsSize( bcNameToBcId.size() );

  computeSourceFluxSizeScalingFactor( time_n,
                                      dt,
                                      domain,
                                      bcNameToBcId,
                                      bcAllSetsSize.toView() );

  // Step 3: we are ready to impose the boundary condition, normalized by the set size

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    integer const isThermal = m_isThermal;

    fsManager.apply< ElementSubRegionBase,
                     SourceFluxBoundaryCondition >( time_n + dt,
                                                    mesh,
                                                    SourceFluxBoundaryCondition::catalogName(),
                                                    [&, isThermal]( SourceFluxBoundaryCondition const & fs,
                                                                    string const & setName,
                                                                    SortedArrayView< localIndex const > const & targetSet,
                                                                    ElementSubRegionBase & subRegion,
                                                                    string const & )
    {
      if( targetSet.size() == 0 )
      {
        return;
      }
      if( !subRegion.hasWrapper( dofKey ) )
      {
        GEOS_LOG_LEVEL_BY_RANK_ON_GROUP( logInfo::SourceFluxFailure,
                                         GEOS_FMT( "{}: trying to apply SourceFlux, but its targetSet named '{}' intersects with non-simulated region named '{}'.",
                                                   getDataContext(), setName, subRegion.getName() ),
                                         fs );
        return;
      }

      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      // Step 3.1: get the values of the source boundary condition that need to be added to the rhs

      array1d< globalIndex > dofArray( targetSet.size() );
      array1d< real64 > rhsContributionArray( targetSet.size() );
      arrayView1d< real64 > rhsContributionArrayView = rhsContributionArray.toView();
      localIndex const rankOffset = dofManager.rankOffset();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > massProd( 0.0 );

      // note that the dofArray will not be used after this step (simpler to use dofNumber instead)
      fs.computeRhsContribution< FieldSpecificationAdd,
                                 parallelDevicePolicy<> >( targetSet.toViewConst(),
                                                           time_n + dt,
                                                           dt,
                                                           subRegion,
                                                           dofNumber,
                                                           rankOffset,
                                                           localMatrix,
                                                           dofArray.toView(),
                                                           rhsContributionArrayView,
                                                           [] GEOS_HOST_DEVICE ( localIndex const )
      {
        return 0.0;
      } );

      // Step 3.2: we are ready to add the right-hand side contributions, taking into account our equation layout

      // get the normalizer
      real64 const sizeScalingFactor = bcAllSetsSize[bcNameToBcId.at( fs.getName())];

      if( isThermal )
      {
        using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
        SingleFluidBase const & fluid =
          getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );

        arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const enthalpy = fluid.enthalpy();
        arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dEnthalpy = fluid.dEnthalpy();
        forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                             targetSet,
                                                             rankOffset,
                                                             ghostRank,
                                                             dofNumber,
                                                             enthalpy,
                                                             dEnthalpy,
                                                             rhsContributionArrayView,
                                                             localRhs,
                                                             localMatrix,
                                                             massProd] GEOS_HOST_DEVICE ( localIndex const a )
        {
          // we need to filter out ghosts here, because targetSet may contain them
          localIndex const ei = targetSet[a];
          if( ghostRank[ei] >= 0 )
          {
            return;
          }

          // add the value to the mass balance equation
          globalIndex const massRowIndex   = dofNumber[ei] - rankOffset;
          globalIndex const energyRowIndex = massRowIndex + 1;
          real64 const rhsValue = rhsContributionArrayView[a] / sizeScalingFactor; // scale the contribution by the sizeScalingFactor here!
          localRhs[massRowIndex] += rhsValue;
          massProd += rhsValue;
          //add the value to the energy balance equation if the flux is positive (i.e., it's a producer)
          if( rhsContributionArrayView[a] > 0.0 )
          {
            globalIndex const pressureDofIndex    = dofNumber[ei] - rankOffset;
            globalIndex const temperatureDofIndex = pressureDofIndex + 1;

            localRhs[energyRowIndex] += enthalpy[ei][0] * rhsValue;

            globalIndex dofIndices[2]{pressureDofIndex, temperatureDofIndex};
            real64 jacobian[2]{rhsValue * dEnthalpy[ei][0][DerivOffset::dP], rhsValue * dEnthalpy[ei][0][DerivOffset::dT]};

            localMatrix.template addToRow< serialAtomic >( energyRowIndex,
                                                           dofIndices,
                                                           jacobian,
                                                           2 );
          }
        } );
      }
      else
      {
        forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                             targetSet,
                                                             rankOffset,
                                                             ghostRank,
                                                             dofNumber,
                                                             rhsContributionArrayView,
                                                             localRhs,
                                                             massProd] GEOS_HOST_DEVICE ( localIndex const a )
        {
          // we need to filter out ghosts here, because targetSet may contain them
          localIndex const ei = targetSet[a];
          if( ghostRank[ei] >= 0 )
          {
            return;
          }

          // add the value to the mass balance equation
          globalIndex const rowIndex = dofNumber[ei] - rankOffset;
          real64 const rhsValue = rhsContributionArrayView[a] / sizeScalingFactor;
          localRhs[rowIndex] += rhsValue;
          massProd += rhsValue;
        } );
      }

      SourceFluxStatsAggregator::forAllFluxStatWrappers( subRegion, fs.getName(),
                                                         [&]( SourceFluxStatsAggregator::WrappedStats & wrapper )
      {
        // set the new sub-region statistics for this timestep
        array1d< real64 > massProdArr{ 1 };
        massProdArr[0] = massProd.get();
        wrapper.gatherTimeStepStats( time_n, dt, massProdArr.toViewConst(), targetSet.size() );
      } );
    } );
  } );
}

void SinglePhaseBase::keepVariablesConstantDuringInitStep( real64 const time,
                                                           real64 const dt,
                                                           DofManager const & dofManager,
                                                           DomainPartition & domain,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time, dt );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      globalIndex const rankOffset = dofManager.rankOffset();
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();

      integer const isThermal = m_isThermal;
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        // 4.1. Apply pressure value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    pres[ei], // freeze the current pressure value
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;

        // 4.2. Apply temperature value to the matrix/rhs
        if( isThermal )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      temp[ei], // freeze the current temperature value
                                                      temp[ei] );
          localRhs[localRow + 1] = rhsValue;
        }
      } );
    } );
  } );
}


void SinglePhaseBase::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );

      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
        updateEnergy( subRegion );
      }
    } );
  } );
}

void SinglePhaseBase::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      arrayView1d< real64 > const pres = subRegion.template getField< fields::flow::pressure >();
      arrayView1d< real64 const > const pres_n = subRegion.template getField< fields::flow::pressure_n >();
      pres.setValues< parallelDevicePolicy<> >( pres_n );

      if( m_isThermal )
      {
        arrayView1d< real64 > const temp = subRegion.template getField< fields::flow::temperature >();
        arrayView1d< real64 const > const temp_n = subRegion.template getField< fields::flow::temperature_n >();
        temp.setValues< parallelDevicePolicy<> >( temp_n );
      }

      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );

      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
        updateEnergy( subRegion );
      }
    } );
  } );
}

real64 SinglePhaseBase::scalingForSystemSolution( DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  real64 scalingFactor = 1.0;
  real64 maxDeltaPres = 0.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      globalIndex const rankOffset = dofManager.rankOffset();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      auto const subRegionData =
        singlePhaseBaseKernels::SolutionScalingKernel::
          launch< parallelDevicePolicy<> >( localSolution, rankOffset, dofNumber, ghostRank, m_maxAbsolutePresChange );

      scalingFactor = std::min( scalingFactor, subRegionData.first );
      maxDeltaPres  = std::max( maxDeltaPres, subRegionData.second );
    } );
  } );

  scalingFactor = MpiWrapper::min( scalingFactor );
  maxDeltaPres  = MpiWrapper::max( maxDeltaPres );

  GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Max pressure change = {} Pa (before scaling)",
                                                      getName(), fmt::format( "{:.{}f}", maxDeltaPres, 3 ) ) );

  return scalingFactor;
}

bool SinglePhaseBase::checkSystemSolution( DomainPartition & domain,
                                           DofManager const & dofManager,
                                           arrayView1d< real64 const > const & localSolution,
                                           real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  integer numNegativePressures = 0;
  real64 minPressure = 0.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      globalIndex const rankOffset = dofManager.rankOffset();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();

      auto const statistics =
        singlePhaseBaseKernels::SolutionCheckKernel::
          launch< parallelDevicePolicy<> >( localSolution, rankOffset, dofNumber, ghostRank, pres, scalingFactor );

      numNegativePressures += statistics.first;
      minPressure = std::min( minPressure, statistics.second );
    } );
  } );

  numNegativePressures = MpiWrapper::sum( numNegativePressures );

  if( numNegativePressures > 0 )
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Number of negative pressure values: {}, minimum value: {} Pa",
                                                        getName(), numNegativePressures, fmt::format( "{:.{}f}", minPressure, 3 ) ) );

  return (m_allowNegativePressure || numNegativePressures == 0) ?  1 : 0;
}

void SinglePhaseBase::saveConvergedState( ElementSubRegionBase & subRegion ) const
{
  FlowSolverBase::saveConvergedState( subRegion );

  arrayView1d< real64 const > const mass = subRegion.template getField< fields::flow::mass >();
  arrayView1d< real64 > const mass_n = subRegion.template getField< fields::flow::mass_n >();
  mass_n.setValues< parallelDevicePolicy<> >( mass );
}

void SinglePhaseBase::applyDeltaVolume( ElementSubRegionBase & subRegion ) const
{
  arrayView1d< real64 > const dVol = subRegion.template getField< fields::flow::deltaVolume >();
  arrayView1d< real64 > const vol = subRegion.template getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString());
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    vol[ei] += dVol[ei];
    dVol[ei] = 0.0;
  } );
}

} /* namespace geos */
