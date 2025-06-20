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
 * @file SinglePhaseBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/AccumulationKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ThermalAccumulationKernels.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidUtils.hpp"


namespace geos
{

namespace constitutive
{

class ConstitutiveBase;

}

/**
 * @class SinglePhaseBase
 *
 * base class to perform a single phase finite volume solve.
 */
class SinglePhaseBase : public FlowSolverBase
{
public:
  using SingleFluidProp = constitutive::SingleFluidVar< real64, 2, constitutive::singlefluid::LAYOUT_FLUID, constitutive::singlefluid::LAYOUT_FLUID_DER >;

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseBase( const string & name,
                   Group * const parent );


  /// deleted default constructor
  SinglePhaseBase() = delete;

  /// deleted copy constructor
  SinglePhaseBase( SinglePhaseBase const & ) = delete;

  /// default move constructor
  SinglePhaseBase( SinglePhaseBase && ) = default;

  /// deleted assignment operator
  SinglePhaseBase & operator=( SinglePhaseBase const & ) = delete;

  /// deleted move operator
  SinglePhaseBase & operator=( SinglePhaseBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseBase() override
  {}

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  ///@{

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  assembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  scalingForSystemSolution( DomainPartition & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  virtual bool
  checkSystemSolution( DomainPartition & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  ///@}

  /**
   * @brief assembles the accumulation terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void assembleAccumulationTerms( DomainPartition & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs );

  /**
   * @brief assembles the accumulation terms for all cells of a spcefici subRegion.
   * @tparam SUBREGION_TYPE the subRegion type
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param subRegion the subRegion
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   *
   */
  template< typename SUBREGION_TYPE >
  void accumulationAssemblyLaunch( DofManager const & dofManager,
                                   SUBREGION_TYPE const & subRegion,
                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                   arrayView1d< real64 > const & localRhs );

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief assembles the flux terms for all cells including jump stabilization
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  virtual void
  assembleStabilizedFluxTerms( real64 const dt,
                               DomainPartition const & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief assembles the flux terms for all cells for the poroelastic case
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   * @param jumpDofKey dofKey of the displacement jump
   */
  virtual void
  assembleEDFMFluxTerms( real64 const time_n,
                         real64 const dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         string const & jumpDofKey ) = 0;

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * elemDofFieldString() { return "singlePhaseVariables"; }
  };

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void
  applyDirichletBC( real64 const time_n,
                    real64 const dt,
                    DomainPartition & domain,
                    DofManager const & dofManager,
                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                    arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply source flux boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void
  applySourceFluxBC( real64 const time_n,
                     real64 const dt,
                     DomainPartition & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply aquifer boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param domain the domain
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  virtual void
  applyAquiferBC( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) const = 0;

  virtual void
  updateState ( DomainPartition & domain ) override final;

  /**
   * @brief Function to update all constitutive state and dependent variables
   * @param subRegion subregion that contains the fields
   */
  real64
  updateFluidState( ElementSubRegionBase & subRegion ) const;


  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  virtual void
  updateFluidModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Function to update fluid mass
   * @param subRegion subregion that contains the fields
   */
  void
  updateMass( ElementSubRegionBase & subRegion ) const;

  /**
   * @brief Template function to update fluid mass
   * @param subRegion subregion that contains the fields
   */
  template< integer IS_THERMAL >
  void updateMass( ElementSubRegionBase & subRegion ) const;

  /**
   * @brief Function to update energy
   * @param subRegion subregion that contains the fields
   */
  void
  updateEnergy( ElementSubRegionBase & subRegion ) const;

  /**
   * @brief Update all relevant solid internal energy models using current values of temperature
   * @param dataGroup the group storing the required fields
   */
  void updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update thermal conductivity
   * @param subRegion the group storing the required fields
   */
  void updateThermalConductivity( ElementSubRegionBase & subRegion ) const;

  /**
   * @brief Function to update fluid mobility
   * @param dataGroup group that contains the fields
   */
  void
  updateMobility( ObjectManagerBase & dataGroup ) const;

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializeFluidState( MeshLevel & mesh, string_array const & regionNames ) override;

  virtual void initializeThermalState( MeshLevel & mesh, string_array const & regionNames ) override;

  /**
   * @brief Compute the hydrostatic equilibrium using the compositions and temperature input tables
   */
  virtual void computeHydrostaticEquilibrium( DomainPartition & domain ) override;

  /**
   * @brief Update the cell-wise pressure gradient
   */
  virtual void updatePressureGradient( DomainPartition & domain )
  {
    GEOS_UNUSED_VAR( domain );
  }

  /**
   * @brief Function to fix the initial state during the initialization step in coupled problems
   * @param[in] time current time
   * @param[in] dt time step
   * @param[in] dofManager degree-of-freedom manager associated with the linear system
   * @param[in] domain the domain
   * @param[in] localMatrix local system matrix
   * @param[in] localRhs local system right-hand side vector
   * @detail This function is meant to be called when the flag m_keepVariablesConstantDuringInitStep is on
   *         The main use case is the initialization step in coupled problems during which we solve an elastic problem for a fixed pressure
   */
  void keepVariablesConstantDuringInitStep( real64 const time,
                                            real64 const dt,
                                            DofManager const & dofManager,
                                            DomainPartition & domain,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs ) const;

  void applyDeltaVolume( ElementSubRegionBase & subRegion ) const;

protected:

  /**
   * @brief Checks constitutive models for consistency
   * @param[in] domain the domain partition
   */
  virtual void validateConstitutiveModels( DomainPartition & domain ) const;

  /**
   * @brief Initialize the aquifer boundary condition (gravity vector, water phase index)
   */
  void initializeAquiferBC() const;

  /**
   * @brief Utility function to save the converged state
   * @param[in] subRegion the element subRegion
   */
  virtual void saveConvergedState( ElementSubRegionBase & subRegion ) const override;

private:

  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

};

template< typename SUBREGION_TYPE >
void SinglePhaseBase::accumulationAssemblyLaunch( DofManager const & dofManager,
                                                  SUBREGION_TYPE const & subRegion,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  if( m_isThermal )
  {
    thermalSinglePhaseBaseKernels::
      AccumulationKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                 dofKey,
                                                 subRegion,
                                                 localMatrix,
                                                 localRhs );
  }
  else
  {
    singlePhaseBaseKernels::
      AccumulationKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                 dofKey,
                                                 subRegion,
                                                 localMatrix,
                                                 localRhs );
  }
}

} /* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_
