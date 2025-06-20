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
 * @file SinglePhasePoromechanics.cpp
 */

#define GEOS_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "SinglePhasePoromechanics.hpp"

#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/solid/PorousDamageSolid.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsDamage.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsLagrangeContact.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVM.hpp"

#include "physicsSolvers/multiphysics/poromechanicsKernels/PoromechanicsKernelsDispatchTypeList.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermoPoromechanicsKernelsDispatchTypeList.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/PoromechanicsDamageKernelsDispatchTypeList.hpp"
#include "physicsSolvers/solidMechanics/kernels/SolidMechanicsKernelsDispatchTypeList.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;
using namespace stabilization;

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::SinglePhasePoromechanics( const string & name,
                                                                                     Group * const parent )
  : Base( name, parent ),
  m_damageFlag()
{
  Base::template addLogLevel< logInfo::LinearSolverConfiguration >();

  this->registerWrapper( viewKeyStruct::damageFlagString(), &m_damageFlag ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The flag to indicate whether a damage solid model is used" );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::postInputInitialization()
{
  Base::postInputInitialization();

  setMGRStrategy();
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::setupCoupling( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                                               DofManager & dofManager ) const
{
  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::setupSystem( DomainPartition & domain,
                                                                             DofManager & dofManager,
                                                                             CRSMatrix< real64, globalIndex > & localMatrix,
                                                                             ParallelVector & rhs,
                                                                             ParallelVector & solution,
                                                                             bool const setSparsity )
{
  if( this->m_precond )
  {
    this->m_precond->clear();
  }

  // setup monolithic coupled system
  PhysicsSolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  dofManager.printFieldInfo();

  if( !this->m_precond && this->m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner();
  }
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::initializePostInitialConditionsPreSubGroups()
{
  Base::initializePostInitialConditionsPreSubGroups();

  string_array const & poromechanicsTargetRegionNames =
    this->template getReference< string_array >( PhysicsSolverBase::viewKeyStruct::targetRegionsString() );
  string_array const & flowTargetRegionNames =
    this->flowSolver()->template getReference< string_array >( PhysicsSolverBase::viewKeyStruct::targetRegionsString() );
  for( size_t i = 0; i < poromechanicsTargetRegionNames.size(); ++i )
  {
    GEOS_THROW_IF( std::find( flowTargetRegionNames.begin(), flowTargetRegionNames.end(), poromechanicsTargetRegionNames[i] )
                   == flowTargetRegionNames.end(),
                   GEOS_FMT( "{} {}: region `{}` must be a target region of `{}`",
                             getCatalogName(), this->getDataContext(), poromechanicsTargetRegionNames[i], this->flowSolver()->getDataContext() ),
                   InputError );
  }
}

template<>
void SinglePhasePoromechanics<>::setMGRStrategy()
{
  LinearSolverParameters & linearSolverParameters = this->m_linearSolverParameters.get();

  if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    return;

  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.dofsPerNode = 3;

  if( dynamic_cast< SinglePhaseHybridFVM * >( this->flowSolver() ) )
  {
    if( this->m_isThermal )
    {
      GEOS_ERROR( GEOS_FMT( "{}: MGR strategy is not implemented for thermal {}/{}",
                            this->getName(), this->getCatalogName(), this->flowSolver()->getCatalogName() ));
    }
    else
    {
      linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics;
    }
  }
  else
  {
    if( this->m_isThermal )
    {
      linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::thermalSinglePhasePoromechanics;
    }
    else
    {
      linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
    }
  }
  GEOS_LOG_LEVEL_RANK_0( logInfo::LinearSolverConfiguration,
                         GEOS_FMT( "{}: MGR strategy set to {}", getName(),
                                   EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy )));
}

template<>
void SinglePhasePoromechanics< SinglePhaseReservoirAndWells<>, SolidMechanicsLagrangianFEM >::setMGRStrategy()
{
  // same as SinglePhaseReservoirAndWells< SinglePhasePoromechanics<> >::setMGRStrategy

  LinearSolverParameters & linearSolverParameters = m_linearSolverParameters.get();

  if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    return;

  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.dofsPerNode = 3;

  if( dynamic_cast< SinglePhaseHybridFVM * >( this->flowSolver() ) )
  {
    GEOS_ERROR( GEOS_FMT( "{}: MGR strategy is not implemented for poromechanics {}/{}",
                          this->getName(), this->getCatalogName(), this->flowSolver()->getCatalogName()));
  }
  else
  {
    linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsReservoirFVM;
  }
  GEOS_LOG_LEVEL_RANK_0( logInfo::LinearSolverConfiguration,
                         GEOS_FMT( "{}: MGR strategy set to {}", this->getName(),
                                   EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy )));
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::assembleSystem( real64 const time_n,
                                                                                real64 const dt,
                                                                                DomainPartition & domain,
                                                                                DofManager const & dofManager,
                                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Steps 1 and 2: compute element-based terms (mechanics and local flow terms)
  assembleElementBasedTerms( time_n,
                             dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  // Step 3: compute the fluxes (face-based contributions)
  if( m_stabilizationType == StabilizationType::Global || m_stabilizationType == StabilizationType::Local )
  {
    this->flowSolver()->assembleStabilizedFluxTerms( dt,
                                                     domain,
                                                     dofManager,
                                                     localMatrix,
                                                     localRhs );
  }
  else
  {
    this->flowSolver()->assembleFluxTerms( dt,
                                           domain,
                                           dofManager,
                                           localMatrix,
                                           localRhs );
  }
}

template<>
void SinglePhasePoromechanics< SinglePhaseReservoirAndWells<>, SolidMechanicsLagrangianFEM >::assembleSystem( real64 const time_n,
                                                                                                              real64 const dt,
                                                                                                              DomainPartition & domain,
                                                                                                              DofManager const & dofManager,
                                                                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                                              arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Steps 1 and 2: compute element-based terms (mechanics and local flow terms)
  assembleElementBasedTerms( time_n,
                             dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  // Step 3: compute the fluxes (face-based contributions)
  if( m_stabilizationType == StabilizationType::Global || m_stabilizationType == StabilizationType::Local )
  {
    this->flowSolver()->assembleStabilizedFluxTerms( dt,
                                                     domain,
                                                     dofManager,
                                                     localMatrix,
                                                     localRhs );
  }
  else
  {
    this->flowSolver()->assembleFluxTerms( dt,
                                           domain,
                                           dofManager,
                                           localMatrix,
                                           localRhs );
  }

  // step 4: assemble well contributions

  this->flowSolver()->wellSolver()->assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
  this->flowSolver()->assembleCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::assembleElementBasedTerms( real64 const time_n,
                                                                                           real64 const dt,
                                                                                           DomainPartition & domain,
                                                                                           DofManager const & dofManager,
                                                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                           arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );

  real64 poromechanicsMaxForce = 0.0;
  real64 mechanicsMaxForce = 0.0;

  // step 1: apply the full poromechanics coupling on the target regions on the poromechanics solver

  set< string > poromechanicsRegionNames;

  this->template forDiscretizationOnMeshTargets<>( domain.getMeshBodies(), [&] ( string const &,
                                                                                 MeshLevel & mesh,
                                                                                 string_array const & regionNames )
  {
    poromechanicsRegionNames.insert( regionNames.begin(), regionNames.end() );

    string const flowDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    if( m_damageFlag )
    {
      poromechanicsMaxForce =
        this->template assemblyLaunch< PoromechanicsDamageKernelsDispatchTypeList,
                                       poromechanicsDamageKernels::SinglePhasePoromechanicsDamageKernelFactory >( mesh,
                                                                                                                  dofManager,
                                                                                                                  regionNames,
                                                                                                                  viewKeyStruct::porousMaterialNamesString(),
                                                                                                                  localMatrix,
                                                                                                                  localRhs,
                                                                                                                  dt,
                                                                                                                  flowDofKey,
                                                                                                                  this->m_performStressInitialization,
                                                                                                                  FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
    else if( this->m_isThermal )
    {
      poromechanicsMaxForce =
        this->template assemblyLaunch< ThermoPoromechanicsKernelsDispatchTypeList,
                                       thermalPoromechanicsKernels::ThermalSinglePhasePoromechanicsKernelFactory >( mesh,
                                                                                                                    dofManager,
                                                                                                                    regionNames,
                                                                                                                    viewKeyStruct::porousMaterialNamesString(),
                                                                                                                    localMatrix,
                                                                                                                    localRhs,
                                                                                                                    dt,
                                                                                                                    flowDofKey,
                                                                                                                    this->m_performStressInitialization,
                                                                                                                    FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
    else
    {
      poromechanicsMaxForce =
        this->template assemblyLaunch< PoromechanicsKernelsDispatchTypeList,
                                       poromechanicsKernels::SinglePhasePoromechanicsKernelFactory >( mesh,
                                                                                                      dofManager,
                                                                                                      regionNames,
                                                                                                      viewKeyStruct::porousMaterialNamesString(),
                                                                                                      localMatrix,
                                                                                                      localRhs,
                                                                                                      dt,
                                                                                                      flowDofKey,
                                                                                                      this->m_performStressInitialization,
                                                                                                      FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
  } );

  // step 2: apply mechanics solver on its target regions not included in the poromechanics solver target regions

  this->solidMechanicsSolver()->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                              MeshLevel & mesh,
                                                                                              string_array const & regionNames )
  {
    // collect the target region of the mechanics solver not included in the poromechanics target regions
    string_array filteredRegionNames;
    filteredRegionNames.reserve( regionNames.size() );
    for( string const & regionName : regionNames )
    {
      // if the mechanics target region is not included in the poromechanics target region, save the string
      if( poromechanicsRegionNames.count( regionName ) == 0 )
      {
        filteredRegionNames.emplace_back( regionName );
      }
    }

    // if the array is empty, the mechanics and poromechanics solver target regions overlap perfectly, there is nothing to do
    if( filteredRegionNames.empty() )
    {
      return;
    }

    mechanicsMaxForce =
      this->template assemblyLaunch< SolidMechanicsKernelsDispatchTypeList,
                                     solidMechanicsLagrangianFEMKernels::QuasiStaticFactory >( mesh,
                                                                                               dofManager,
                                                                                               filteredRegionNames,
                                                                                               SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                                                               localMatrix,
                                                                                               localRhs,
                                                                                               dt );
  } );

  this->solidMechanicsSolver()->applyContactConstraint( dofManager, domain, localMatrix, localRhs );
  this->solidMechanicsSolver()->getMaxForce() = LvArray::math::max( mechanicsMaxForce, poromechanicsMaxForce );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::createPreconditioner()
{
  if( this->m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( this->solidMechanicsSolver()->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { solidMechanics::totalDisplacement::key(), { 3, true } } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( this->flowSolver()->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { flow::pressure::key(), { 1, true } } },
                         std::move( flowPrecond ) );

    this->m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::updateBulkDensity( ElementSubRegionBase & subRegion )
{
  // get the fluid model (to access fluid density)
  string const fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
  SingleFluidBase const & fluid = this->template getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

  // get the solid model (to access porosity and solid density)
  string const solidName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
  CoupledSolidBase const & solid = this->template getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

  // update the bulk density
  poromechanicsKernels::
    SinglePhaseBulkDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( fluid,
                                               solid,
                                               subRegion );
}

template class SinglePhasePoromechanics<>;
template class SinglePhasePoromechanics< SinglePhaseBase, SolidMechanicsLagrangeContact >;
template class SinglePhasePoromechanics< SinglePhaseBase, SolidMechanicsEmbeddedFractures >;
template class SinglePhasePoromechanics< SinglePhaseReservoirAndWells<> >;
template class SinglePhasePoromechanics< SinglePhaseReservoirAndWells<>, SolidMechanicsLagrangeContact >;
//template class SinglePhasePoromechanics< SinglePhaseReservoirAndWells<>, SolidMechanicsEmbeddedFractures >;

namespace
{
typedef SinglePhasePoromechanics< SinglePhaseReservoirAndWells<> > SinglePhaseReservoirPoromechanics;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseReservoirPoromechanics, string const &, Group * const )
typedef SinglePhasePoromechanics<> SinglePhasePoromechanics;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhasePoromechanics, string const &, Group * const )
}

} /* namespace geos */
