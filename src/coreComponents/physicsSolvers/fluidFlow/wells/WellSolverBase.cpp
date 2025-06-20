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
 * @file WellSolverBase.cpp
 */

#include "WellSolverBase.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/ThermalCompositionalMultiphaseWellKernels.hpp"
#include "fileIO/Outputs/OutputBase.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

WellSolverBase::WellSolverBase( string const & name,
                                Group * const parent )
  : PhysicsSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_numDofPerWellElement( 0 ),
  m_numDofPerResElement( 0 ),
  m_isThermal( 0 ),
  m_ratesOutputDir( joinPath( OutputBase::getOutputDirectory(), name + "_rates" ) ),
  m_keepVariablesConstantDuringInitStep( false )
{
  registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not." );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

  this->registerWrapper( viewKeyStruct::writeCSVFlagString(), &m_writeCSV ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "When set to 1, write the rates into a CSV file" );

  this->registerWrapper( viewKeyStruct::timeStepFromTablesFlagString(), &m_timeStepFromTables ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Choose time step to honor rates/bhp tables time intervals" );
}

Group *WellSolverBase::createChild( string const & childKey, string const & childName )
{
  const auto childTypes = { keys::wellControls };
  GEOS_ERROR_IF( childKey != keys::wellControls,
                 PhysicsSolverBase::CatalogInterface::unknownTypeError( childKey, getDataContext(), childTypes ) );
  return &registerGroup< WellControls >( childName );
}

void WellSolverBase::expandObjectCatalogs()
{
  createChild( keys::wellControls, keys::wellControls );
}

WellSolverBase::~WellSolverBase() = default;

void WellSolverBase::postInputInitialization()
{
  PhysicsSolverBase::postInputInitialization();

  // 1. Set key dimensions of the problem
  m_numDofPerWellElement = m_isThermal ?    m_numComponents + 2 : m_numComponents + 1; // 1 pressure  connectionRate + temp if thermal
  m_numDofPerResElement = m_isThermal ? m_numComponents  + 1: m_numComponents;   // 1 pressure   + temp if thermal


  // create dir for rates output
  if( m_writeCSV > 0 )
  {
    if( MpiWrapper::commRank() == 0 )
    {
      makeDirsForPath( m_ratesOutputDir );
    }
    // wait till the dir is created by rank 0
    MPI_Barrier( MPI_COMM_WORLD );
  }
}

void WellSolverBase::registerDataOnMesh( Group & meshBodies )
{
  PhysicsSolverBase::registerDataOnMesh( meshBodies );

  // loop over the wells
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    string_array const & regionNames )
  {

    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            WellElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::well::pressure >( getName() );
      subRegion.registerField< fields::well::pressure_n >( getName() );

      subRegion.registerField< fields::well::temperature >( getName() );
      if( isThermal() )
      {
        subRegion.registerField< fields::well::temperature_n >( getName() );
      }

      subRegion.registerField< fields::well::gravityCoefficient >( getName() );

      subRegion.registerWrapper< string >( viewKeyStruct::fluidNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );

      PerforationData * const perforationData = subRegion.getPerforationData();
      perforationData->registerField< fields::well::gravityCoefficient >( getName() );
    } );
  } );
}

void WellSolverBase::initializePostSubGroups()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      validateWellConstraints( 0, 0, subRegion );
    } );
  } );
}

void WellSolverBase::setupDofs( DomainPartition const & domain,
                                DofManager & dofManager ) const
{
  map< std::pair< string, string >, string_array > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                               MeshLevel const & meshLevel,
                                                               string_array const & regionNames )
  {
    string_array regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      WellElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    auto const key = std::make_pair( meshBodyName, meshLevel.getName());
    meshTargets[key] = std::move( regions );
  } );

  dofManager.addField( wellElementDofName(),
                       FieldLocation::Elem,
                       numDofPerWellElement(),
                       meshTargets );

  dofManager.addCoupling( wellElementDofName(),
                          wellElementDofName(),
                          DofManager::Connector::Node );
}

void WellSolverBase::implicitStepSetup( real64 const & time_n,
                                        real64 const & GEOS_UNUSED_PARAM( dt ),
                                        DomainPartition & domain )
{
  // Initialize the primary and secondary variables for the first time step

  initializeWells( domain, time_n );
}

void WellSolverBase::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    { updateSubRegionState( subRegion ); } );
  } );
}

void WellSolverBase::assembleSystem( real64 const time,
                                     real64 const dt,
                                     DomainPartition & domain,
                                     DofManager const & dofManager,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs )
{
  string const wellDofKey = dofManager.getKey( wellElementDofName());

  // assemble the accumulation term in the mass balance equations
  assembleAccumulationTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // then assemble the pressure relations between well elements
  assemblePressureRelations( time, dt, domain, dofManager, localMatrix, localRhs );
  // then compute the perforation rates (later assembled by the coupled solver)
  computePerforationRates( time, dt, domain );

  // then assemble the flux terms in the mass balance equations
  // get a reference to the degree-of-freedom numbers
  // then assemble the flux terms in the mass balance equations
  assembleFluxTerms( time, dt, domain, dofManager, localMatrix, localRhs );
}

void WellSolverBase::initializePostInitialConditionsPreSubGroups()
{
  PhysicsSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // make sure that nextWellElementIndex is up-to-date (will be used in well initialization and assembly)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    { subRegion.reconstructLocalConnectivity(); } );
  } );

  // Precompute solver-specific constant data (e.g. gravity-coefficient)
  precomputeData( domain );
}

void WellSolverBase::precomputeData( DomainPartition & domain )
{
  R1Tensor const gravVector = gravityVector();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      PerforationData & perforationData = *subRegion.getPerforationData();
      WellControls & wellControls = getWellControls( subRegion );
      real64 const refElev = wellControls.getReferenceElevation();

      arrayView2d< real64 const > const wellElemLocation = subRegion.getElementCenter();
      arrayView1d< real64 > const wellElemGravCoef = subRegion.getField< fields::well::gravityCoefficient >();

      arrayView2d< real64 const > const perfLocation = perforationData.getField< fields::perforation::location >();
      arrayView1d< real64 > const perfGravCoef = perforationData.getField< fields::well::gravityCoefficient >();

      forAll< serialPolicy >( perforationData.size(), [=]( localIndex const iperf )
      {
        // precompute the depth of the perforations
        perfGravCoef[iperf] = LvArray::tensorOps::AiBi< 3 >( perfLocation[iperf], gravVector );
      } );

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
      {
        // precompute the depth of the well elements
        wellElemGravCoef[iwelem] = LvArray::tensorOps::AiBi< 3 >( wellElemLocation[iwelem], gravVector );
      } );

      // set the reference well element where the BHP control is applied
      wellControls.setReferenceGravityCoef( refElev * gravVector[2] );
    } );
  } );
}

WellControls & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion )
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

WellControls const & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion ) const
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

real64 WellSolverBase::setNextDt( real64 const & currentTime, const real64 & currentDt, geos::DomainPartition & domain )
{
  real64 nextDt = PhysicsSolverBase::setNextDt( currentTime, currentDt, domain );

  if( m_timeStepFromTables )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            WellElementSubRegion & subRegion )
      {
        WellControls & wellControls = getWellControls( subRegion );
        real64 const nextDt_orig = nextDt;
        wellControls.setNextDtFromTables( currentTime, nextDt );
        if( m_nonlinearSolverParameters.getLogLevel() > 0 && nextDt < nextDt_orig )
          GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on tables coordinates = {}", getName(), nextDt ));
      } );
    } );
  }

  return nextDt;
}

} // namespace geos
