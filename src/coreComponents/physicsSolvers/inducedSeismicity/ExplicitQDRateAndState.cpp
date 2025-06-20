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
 * @file ExplicitQDRateAndState.cpp
 */

#include "ExplicitQDRateAndState.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mesh/DomainPartition.hpp"
#include "kernels/ExplicitRateAndStateKernels.hpp"
#include "kernels/EmbeddedRungeKuttaKernels.hpp"
#include "rateAndStateFields.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;
using namespace rateAndStateKernels;

ExplicitQDRateAndState::ExplicitQDRateAndState( const string & name,
                                                Group * const parent ):
  QDRateAndStateBase( name, parent ),
  m_butcherTable( BogackiShampine32Table()), // TODO: The butcher table should be specified in the XML input.
  m_successfulStep( false ),
  m_stepUpdateFactor( 1.0 ),
  m_controller( PIDController( { 0.6, -0.2, 0.0 },
                               1.0e-6, 1.0e-6, 0.81 )) // TODO: The control parameters should be specified in the XML input
{
  addLogLevel< logInfo::SolverSteps >();
}

ExplicitQDRateAndState::~ExplicitQDRateAndState()
{
  // TODO Auto-generated destructor stub
}

void ExplicitQDRateAndState::registerDataOnMesh( Group & meshBodies )
{
  QDRateAndStateBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      SurfaceElementSubRegion & subRegion )
    {
      // Runge-Kutta stage rates and error
      integer const numRKComponents = 3;
      subRegion.registerField< rateAndState::rungeKuttaStageRates >( getName() ).reference().resizeDimension< 1, 2 >( m_butcherTable.numStages, numRKComponents );
      subRegion.registerField< rateAndState::error >( getName() ).reference().resizeDimension< 1 >( numRKComponents );
    } );
  } );
}

real64 ExplicitQDRateAndState::solverStep( real64 const & time_n,
                                           real64 const & dt,
                                           int const cycleNumber,
                                           DomainPartition & domain )
{
  applyInitialConditionsToFault( cycleNumber, domain );
  saveState( domain );

  real64 dtAdaptive = dt;

  GEOS_LOG_LEVEL_RANK_0( logInfo::SolverSteps, "Rate and State solver" );
  while( true ) // Adaptive time step loop. Performs a Runge-Kutta time stepping with error control on state and slip
  {
    real64 dtStress; GEOS_UNUSED_VAR( dtStress );

    // Initial Runge-Kutta stage
    stepRateStateODEInitialSubstage( dtAdaptive, domain );
    real64 dtStage = m_butcherTable.c[1]*dtAdaptive;
    dtStress = updateStresses( time_n, dtStage, cycleNumber, domain );
    updateSlipVelocity( time_n, dtStage, domain );

    // Remaining stages
    for( integer stageIndex = 1; stageIndex < m_butcherTable.numStages-1; stageIndex++ )
    {
      stepRateStateODESubstage( stageIndex, dtAdaptive, domain );
      dtStage = m_butcherTable.c[stageIndex+1]*dtAdaptive;
      dtStress = updateStresses( time_n, dtStage, cycleNumber, domain );
      updateSlipVelocity( time_n, dtStage, domain );
    }

    stepRateStateODEAndComputeError( dtAdaptive, domain );
    // Update timestep based on the time step error
    evalTimestep( domain );
    if( m_successfulStep ) // set in evalTimestep
    {
      // Compute stresses, and slip velocity and save results at time_n + dtAdapitve
      if( !m_butcherTable.FSAL )
      {
        dtStress = updateStresses( time_n, dtAdaptive, cycleNumber, domain );
        updateSlipVelocity( time_n, dtAdaptive, domain );
      }
      saveState( domain );
      break;
    }
    else
    {
      // Retry with updated time step
      dtAdaptive = setNextDt( time_n, dtAdaptive, domain );
    }
  }
  // return last successful adaptive time step (passed along to setNextDt)
  return dtAdaptive;
}

void ExplicitQDRateAndState::stepRateStateODEInitialSubstage( real64 const dt, DomainPartition & domain ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {

      string const & frictionLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
      constitutive::ConstitutiveBase & frictionLaw = subRegion.getConstitutiveModel< constitutive::ConstitutiveBase >( frictionLawName );
      constitutive::ConstitutivePassThru< constitutive::RateAndStateFrictionBase >::execute( frictionLaw, [&] ( auto & castedFrictionLaw )
      {
        rateAndStateKernels::createAndlaunchODEInitialSubStage( subRegion, castedFrictionLaw, m_butcherTable, dt, m_successfulStep );
      } );
    } );
  } );
}

void ExplicitQDRateAndState::stepRateStateODESubstage( integer const stageIndex,
                                                       real64 const dt,
                                                       DomainPartition & domain ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {

      string const & frictionLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
      constitutive::ConstitutiveBase & frictionLaw = subRegion.getConstitutiveModel< constitutive::ConstitutiveBase >( frictionLawName );
      constitutive::ConstitutivePassThru< constitutive::RateAndStateFrictionBase >::execute( frictionLaw, [&] ( auto & castedFrictionLaw )
      {
        rateAndStateKernels::createAndlaunchStepRateStateODESubstage( subRegion, castedFrictionLaw, m_butcherTable, stageIndex, dt );
      } );
    } );
  } );
}

void ExplicitQDRateAndState::stepRateStateODEAndComputeError( real64 const dt, DomainPartition & domain ) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {

      string const & frictionLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
      constitutive::ConstitutiveBase & frictionLaw = getConstitutiveModel< constitutive::ConstitutiveBase >( subRegion, frictionLawName );
      constitutive::ConstitutivePassThru< constitutive::RateAndStateFrictionBase >::execute( frictionLaw, [&] ( auto & castedFrictionLaw )
      {
        rateAndStateKernels::createAndlaunchStepRateStateODEAndComputeError( subRegion,
                                                                             castedFrictionLaw,
                                                                             m_butcherTable,
                                                                             m_controller.relTol,
                                                                             m_controller.absTol,
                                                                             dt );
      } );
    } );
  } );
}

void ExplicitQDRateAndState::updateSlipVelocity( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain ) const
{
  GEOS_LOG_LEVEL_RANK_0( logInfo::SolverSteps, "Rate and State solver" );
  integer const maxIterNewton = m_nonlinearSolverParameters.m_maxIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {
      string const & frictionLawName = subRegion.getReference< string >( viewKeyStruct::frictionLawNameString() );
      constitutive::ConstitutiveBase & frictionLaw = subRegion.getConstitutiveModel< constitutive::ConstitutiveBase >( frictionLawName );
      constitutive::ConstitutivePassThru< constitutive::RateAndStateFrictionBase >::execute( frictionLaw, [=, &subRegion] ( auto & castedFrictionLaw )
      {
        // solve rate and state equations.
        rateAndStateKernels::createAndLaunch< rateAndStateKernels::ExplicitRateAndStateKernel,
                                              parallelDevicePolicy<> >( subRegion,
                                                                        castedFrictionLaw,
                                                                        m_shearImpedance,
                                                                        maxIterNewton,
                                                                        newtonTol,
                                                                        time_n,
                                                                        dt );
      } );
    } );
  } );
}

void ExplicitQDRateAndState::evalTimestep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion const & subRegion )
    {
      arrayView2d< real64 const > const error = subRegion.getField< rateAndState::error >();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > scaledl2ErrorSquared( 0.0 );
      integer const N = subRegion.size();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        scaledl2ErrorSquared += LvArray::tensorOps::l2NormSquared< 3 >( error[k] );
      } );
      m_controller.errors[0] = LvArray::math::sqrt( MpiWrapper::sum( scaledl2ErrorSquared.get() / (3.0*N) ));
    } );
  } );

  // Compute update factor to currentDt using PID error controller + limiter
  m_stepUpdateFactor = m_controller.computeUpdateFactor( m_butcherTable.algHighOrder, m_butcherTable.algLowOrder );
  // Check if step was acceptable
  m_successfulStep = (m_stepUpdateFactor >= m_controller.acceptSafety) ? true : false;
  if( m_successfulStep )
  {
    m_controller.errors[2] = m_controller.errors[1];
    m_controller.errors[1] = m_controller.errors[0];
  }
}

real64 ExplicitQDRateAndState::setNextDt( real64 const & currentTime,
                                          real64 const & currentDt,
                                          DomainPartition & domain )
{
  GEOS_UNUSED_VAR( currentTime, domain );
  real64 const nextDt = m_stepUpdateFactor*currentDt;
  if( m_successfulStep )
  {
    GEOS_LOG_LEVEL_RANK_0( logInfo::SolverSteps,
                           GEOS_FMT( "Adaptive time step successful. The next dt will be {:.2e} s", nextDt ));
  }
  else
  {
    GEOS_LOG_LEVEL_RANK_0( logInfo::SolverSteps,
                           GEOS_FMT( "Adaptive time step failed. Retry step with dt {:.2e} s", nextDt ));
  }
  return nextDt;
}

} // namespace geos
