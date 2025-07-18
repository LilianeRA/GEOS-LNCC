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
 * @file SinglePhasePoromechanicsEmbeddedFractures.cpp
 */

#include "SinglePhasePoromechanicsEmbeddedFractures.hpp"
#include "constitutive/contact/HydraulicApertureRelationSelector.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsEFEMKernelsHelper.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEFEM.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsEFEM.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsFractures.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"


namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;

SinglePhasePoromechanicsEmbeddedFractures::SinglePhasePoromechanicsEmbeddedFractures( const std::string & name,
                                                                                      Group * const parent ):
  SinglePhasePoromechanics( name, parent )
{
  Base::template addLogLevel< logInfo::LinearSolverConfiguration >();
}

SinglePhasePoromechanicsEmbeddedFractures::~SinglePhasePoromechanicsEmbeddedFractures()
{}

void SinglePhasePoromechanicsEmbeddedFractures::setMGRStrategy()
{
  LinearSolverParameters & linearSolverParameters = m_linearSolverParameters.get();

  if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    return;

  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.dofsPerNode = 3;

  linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsEmbeddedFractures;
  GEOS_LOG_LEVEL_RANK_0( logInfo::LinearSolverConfiguration,
                         GEOS_FMT( "{}: MGR strategy set to {}", getName(),
                                   EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy )));
}

void SinglePhasePoromechanicsEmbeddedFractures::postInputInitialization()
{
  Base::postInputInitialization();

  GEOS_ERROR_IF( solidMechanicsSolver()->useStaticCondensation(),
                 GEOS_FMT( "{}: {} = 1 in {} solver named {} is not supported for {}",
                           this->getName(), SolidMechanicsEmbeddedFractures::viewKeyStruct::useStaticCondensationString(),
                           solidMechanicsSolver()->getCatalogName(), solidMechanicsSolver()->getName(), getCatalogName() ));
}

void SinglePhasePoromechanicsEmbeddedFractures::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  Base::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&] ( localIndex const,
                                                                                     EmbeddedSurfaceSubRegion & subRegion )
    {
      subRegion.registerField< contact::dTraction_dPressure >( getName() );
    } );
  } );
}

void SinglePhasePoromechanicsEmbeddedFractures::initializePostInitialConditionsPreSubGroups()
{
  Base::initializePostInitialConditionsPreSubGroups();

  updateState( this->getGroupByPath< DomainPartition >( "/Problem/domain" ) );
}

void SinglePhasePoromechanicsEmbeddedFractures::setupCoupling( DomainPartition const & domain,
                                                               DofManager & dofManager ) const
{
  Base::setupCoupling( domain, dofManager );

  map< std::pair< string, string >, string_array > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel const & meshLevel,
                                                                string_array const & regionNames )
  {
    string_array regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< SurfaceElementRegion >( regionNames,
                                                                    [&]( localIndex const,
                                                                         SurfaceElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    meshTargets[std::make_pair( meshBodyName, meshLevel.getName())] = std::move( regions );
  } );

  dofManager.addCoupling( SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          contact::dispJump::key(),
                          DofManager::Connector::Elem,
                          meshTargets );
}

void SinglePhasePoromechanicsEmbeddedFractures::setupSystem( DomainPartition & domain,
                                                             DofManager & dofManager,
                                                             CRSMatrix< real64, globalIndex > & localMatrix,
                                                             ParallelVector & rhs,
                                                             ParallelVector & solution,
                                                             bool const setSparsity )
{
  // Add missing couplings ( matrix pressure with displacement jump and jump - displacement )

  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( setSparsity );

  dofManager.setDomain( domain );
  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  // Set the sparsity pattern without the Kwu and Kuw blocks.
  SparsityPattern< globalIndex > patternDiag;
  dofManager.setSparsityPattern( patternDiag );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternDiag.numRows() );
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternDiag.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling jump-pm
  addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ) );
  }

  dofManager.printFieldInfo();

  // Add the nonzeros from coupling
  addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.setName( this->getName() + "/localMatrix" );
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );
}

void SinglePhasePoromechanicsEmbeddedFractures::addCouplingNumNonzeros( DomainPartition & domain,
                                                                        DofManager & dofManager,
                                                                        arrayView1d< localIndex > const & rowLengths ) const
{
  // 1. Add the number of nonzeros induced by coupling jump-displacement
  solidMechanicsSolver()->addCouplingNumNonzeros( domain, dofManager, rowLengths );

  // 2. Add the number of nonzeros induced by coupling jump - matrix flow dofs
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )

  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const jumpDofKey = dofManager.getKey( contact::dispJump::key() );
    string const flowDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    globalIndex const rankOffset = dofManager.rankOffset();

    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion )
    {
      localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();

      OrderedVariableToManyElementRelation const & embeddedSurfacesToCells = embeddedSurfaceSubRegion.getToCellRelation();

      arrayView1d< globalIndex const > const &
      embeddedElementDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );
      arrayView1d< integer const > const & ghostRank = embeddedSurfaceSubRegion.ghostRank();

      for( localIndex k=0; k<numEmbeddedElems; ++k )
      {
        // Get rock matrix element subregion
        CellElementSubRegion const & subRegion =
          elemManager.getRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] ).
            getSubRegion< CellElementSubRegion >( embeddedSurfacesToCells.m_toElementSubRegion[k][0] );

        arrayView1d< globalIndex const > const &
        flowDofNumber = subRegion.getReference< globalIndex_array >( flowDofKey );

        localIndex cellElementIndex = embeddedSurfacesToCells.m_toElementIndex[k][0];

        if( ghostRank[k] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( embeddedElementDofNumber[k] - rankOffset );
          GEOS_ASSERT_GE( localRow, 0 );
          GEOS_ASSERT_GE( rowLengths.size(), localRow + embeddedSurfaceSubRegion.numOfJumpEnrichments()  );

          for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); ++i )
          {
            rowLengths[localRow + i] += flowSolver()->numberOfDofsPerCell();
          }
          for( integer i = 0; i < flowSolver()->numberOfDofsPerCell(); i++ )
          {
            localIndex const localFlowDofRow = LvArray::integerConversion< localIndex >( flowDofNumber[cellElementIndex] + i - rankOffset );
            GEOS_ASSERT_GE( localFlowDofRow, 0 );
            GEOS_ASSERT_GE( rowLengths.size(), localFlowDofRow + embeddedSurfaceSubRegion.numOfJumpEnrichments() );

            rowLengths[ localFlowDofRow ] += embeddedSurfaceSubRegion.numOfJumpEnrichments();
          }
        }
      }
    } );

    // 3. Add the number of nonzeros induced by coupling jump (aperture) - fracture pressure due to flux term
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( flowSolver()->getDiscretizationName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
          elemManager.getRegion( seri[iconn][0] ).getSubRegion< EmbeddedSurfaceSubRegion >( sesri[iconn][0] );

        arrayView1d< globalIndex const > const &
        flowDofNumber =  embeddedSurfaceSubRegion.getReference< globalIndex_array >( flowDofKey );

        for( localIndex k0=0; k0<numFluxElems; ++k0 )
        {
          globalIndex const activeFlowDOF = flowDofNumber[sei[iconn][k0]];
          globalIndex const rowNumber = activeFlowDOF - rankOffset;

          if( rowNumber >= 0 && rowNumber < rowLengths.size() )
          {
            for( localIndex k1=0; k1<numFluxElems; ++k1 )
            {
              // The coupling with the jump of the cell itself has already been added by the dofManager
              // so we only add the coupling with the jumps of the neighbours.
              if( k1 != k0 )
              {
                rowLengths[ rowNumber ] += embeddedSurfaceSubRegion.numOfJumpEnrichments(); // number of jump enrichments.
                if( m_isThermal )
                {
                  // energy flux is also coupled to dispJump
                  rowLengths[ rowNumber + 1 ] += embeddedSurfaceSubRegion.numOfJumpEnrichments();
                }
              }
            }
          }
        }
      }
    } );
  } );
}

void SinglePhasePoromechanicsEmbeddedFractures::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                            DofManager const & dofManager,
                                                                            SparsityPatternView< globalIndex > const & pattern ) const
{
  // 1. Add sparsity pattern induced by coupling jump-displacement
  solidMechanicsSolver()->addCouplingSparsityPattern( domain, dofManager, pattern );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const jumpDofKey = dofManager.getKey( fields::contact::dispJump::key() );
    string const flowDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    globalIndex const rankOffset = dofManager.rankOffset();

    // 2. Add the sparsity pattern induced by coupling jump - matrix pressure
    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&]( localIndex const, EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion )
    {
      localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();

      OrderedVariableToManyElementRelation const & embeddedSurfacesToCells = embeddedSurfaceSubRegion.getToCellRelation();

      arrayView1d< globalIndex const > const &
      jumpDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );
      arrayView1d< integer const > const & ghostRank = embeddedSurfaceSubRegion.ghostRank();

      for( localIndex k=0; k<numEmbeddedElems; ++k )
      {
        // Get rock matrix element subregion
        CellElementSubRegion const & subRegion =
          elemManager.getRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] ).
            getSubRegion< CellElementSubRegion >( embeddedSurfacesToCells.m_toElementSubRegion[k][0] );

        arrayView1d< globalIndex const > const &
        flowDofNumber = subRegion.getReference< globalIndex_array >( flowDofKey );

        localIndex cellElementIndex = embeddedSurfacesToCells.m_toElementIndex[k][0];

        if( ghostRank[k] < 0 ) /// TODO is this really necessary?
        {
          localIndex const localJumpRow = LvArray::integerConversion< localIndex >( jumpDofNumber[k] - rankOffset );
          for( integer dof = 0; dof  < flowSolver()->numberOfDofsPerCell(); dof++ )
          {

            localIndex const localFlowDofRow = LvArray::integerConversion< localIndex >( flowDofNumber[cellElementIndex] + dof - rankOffset );

            for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); ++i )
            {
              if( localJumpRow + i >= 0 && localJumpRow + i < pattern.numRows() )
                pattern.insertNonZero( localJumpRow + i, flowDofNumber[cellElementIndex] + dof );
              if( localFlowDofRow >= 0 && localFlowDofRow < pattern.numRows() )
                pattern.insertNonZero( localFlowDofRow, jumpDofNumber[k] + i );
            }
          }
        }
      }
    } );

    // 3. Add the sparsity pattern induced by coupling jump (aperture) - fracture pressure due to flux term
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( flowSolver()->getDiscretizationName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
          elemManager.getRegion( seri[iconn][0] ).getSubRegion< EmbeddedSurfaceSubRegion >( sesri[iconn][0] );

        arrayView1d< globalIndex const > const &
        flowDofNumber =  embeddedSurfaceSubRegion.getReference< globalIndex_array >( flowDofKey );
        arrayView1d< globalIndex const > const &
        jumpDofNumber =  embeddedSurfaceSubRegion.getReference< globalIndex_array >( jumpDofKey );

        for( localIndex k0=0; k0<numFluxElems; ++k0 )
        {
          for( integer dof = 0; dof < flowSolver()->numberOfDofsPerCell(); dof++ )
          {
            globalIndex const activeFlowDOF = flowDofNumber[sei[iconn][k0]] + dof;
            localIndex const rowIndex = activeFlowDOF - rankOffset;

            if( rowIndex >= 0 && rowIndex < pattern.numRows() )
            {
              for( localIndex k1=0; k1<numFluxElems; ++k1 )
              {
                // The coupling with the jump of the cell itself has already been added by the dofManager
                // so we only add the coupling with the jumps of the neighbours.
                if( k1 != k0 )
                {
                  for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); i++ )
                  {
                    globalIndex const colIndex = jumpDofNumber[sei[iconn][k1]] + i;
                    pattern.insertNonZero( rowIndex, colIndex );
                  }
                }
              }
            }
          }
        }
      }
    } );
  } );

}

void SinglePhasePoromechanicsEmbeddedFractures::assembleSystem( real64 const time_n,
                                                                real64 const dt,
                                                                DomainPartition & domain,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  //updateState( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )

  {
    if( m_isThermal )
    {
      solidMechanicsSolver()->getMaxForce() =
        assemblyLaunch< thermalPoromechanicsKernels::ThermalSinglePhasePoromechanicsKernelFactory,
                        thermoPoromechanicsEFEMKernels::ThermalSinglePhasePoromechanicsEFEMKernelFactory >( mesh,
                                                                                                            dofManager,
                                                                                                            regionNames,
                                                                                                            Base::viewKeyStruct::porousMaterialNamesString(),
                                                                                                            localMatrix,
                                                                                                            localRhs,
                                                                                                            dt );
    }
    else
    {
      solidMechanicsSolver()->getMaxForce() =
        assemblyLaunch< poromechanicsKernels::SinglePhasePoromechanicsKernelFactory,
                        poromechanicsEFEMKernels::SinglePhaseKernelFactory >( mesh,
                                                                              dofManager,
                                                                              regionNames,
                                                                              Base::viewKeyStruct::porousMaterialNamesString(),
                                                                              localMatrix,
                                                                              localRhs,
                                                                              dt );
    }

    // 3. Assemble poroelastic fluxes and all derivatives
    string const jumpDofKey = dofManager.getKey( contact::dispJump::key() );
    flowSolver()->assembleEDFMFluxTerms( time_n, dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs,
                                         jumpDofKey );

  } );

}

void SinglePhasePoromechanicsEmbeddedFractures::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  /// 1. update the reservoir
  Base::updateState( domain );

  // remove the contribution of the hydraulic aperture from the stencil weights
  flowSolver()->prepareStencilWeights( domain );

  /// 2. update the fractures
  solidMechanicsSolver()->updateState( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&] ( localIndex const,
                                                                                     auto & subRegion )
    {
      arrayView2d< real64 const > const dispJump =
        subRegion.template getField< contact::dispJump >();

      arrayView1d< real64 > const aperture = subRegion.getElementAperture();

      arrayView1d< real64 > const hydraulicAperture =
        subRegion.template getField< flow::hydraulicAperture >();

      arrayView1d< real64 const > const oldHydraulicAperture =
        subRegion.template getField< flow::aperture0 >();

      arrayView1d< real64 const > const volume = subRegion.getElementVolume();

      arrayView1d< real64 > const deltaVolume =
        subRegion.template getField< flow::deltaVolume >();

      arrayView1d< real64 const > const area = subRegion.getElementArea().toViewConst();

      arrayView2d< real64 > const & fractureContactTraction = subRegion.template getField< contact::traction >();

      arrayView1d< real64 const > const & pressure =
        subRegion.template getField< flow::pressure >();

      string const & hydraulicApertureRelationName = subRegion.template getReference< string >( viewKeyStruct::hydraulicApertureRelationNameString()  );
      HydraulicApertureBase const & hydraulicApertureModel = this->template getConstitutiveModel< HydraulicApertureBase >( subRegion, hydraulicApertureRelationName );

      string const porousSolidName = subRegion.template getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
      CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( porousSolidName );

      constitutive::ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion, &hydraulicApertureModel] ( auto & castedPorousSolid )
      {
        typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousMaterialWrapper = castedPorousSolid.createKernelUpdates();

        constitutiveUpdatePassThru( hydraulicApertureModel, [=, &subRegion] ( auto & castedHydraulicApertureModel )
        {

          using HydraulicApertureModelType = TYPEOFREF( castedHydraulicApertureModel );
          typename HydraulicApertureModelType::KernelWrapper hydraulicApertureModelWrapper = castedHydraulicApertureModel.createKernelWrapper();

          poromechanicsFracturesKernels::StateUpdateKernel::
            launch< parallelDevicePolicy<> >( subRegion.size(),
                                              porousMaterialWrapper,
                                              hydraulicApertureModelWrapper,
                                              dispJump,
                                              pressure,
                                              area,
                                              volume,
                                              deltaVolume,
                                              aperture,
                                              oldHydraulicAperture,
                                              hydraulicAperture,
                                              fractureContactTraction );

        } );
      } );

      // update the stencil weights using the updated hydraulic aperture
      flowSolver()->updateStencilWeights( domain );
      // update fracture's porosity from pressure and temperature
      flowSolver()->updatePorosityAndPermeability( subRegion );
      // update fluid model
      flowSolver()->updateFluidState( subRegion );
      if( m_isThermal )
      {
        // update solid internal energy
        flowSolver()->updateSolidInternalEnergyModel( subRegion );
      }
    } );
  } );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhasePoromechanicsEmbeddedFractures, std::string const &, Group * const )

} /* namespace geos */
