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
 * @file SinglePhasePoromechanicsEFEM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEFEM_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEFEM_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace poromechanicsEFEMKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhasePoromechanicsEFEM :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 0 >;
  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  /// Compile time value for the number of gotquadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;


  SinglePhasePoromechanicsEFEM( NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager const & faceManager,
                                localIndex const targetRegionIndex,
                                SUBREGION_TYPE const & elementSubRegion,
                                FE_TYPE const & finiteElementSpace,
                                CONSTITUTIVE_TYPE & inputConstitutiveType,
                                EmbeddedSurfaceSubRegion const & embeddedSurfSubRegion,
                                arrayView1d< globalIndex const > const dispDofNumber,
                                arrayView1d< globalIndex const > const jumpDofNumber,
                                string const inputFlowDofKey,
                                globalIndex const rankOffset,
                                CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                arrayView1d< real64 > const inputRhs,
                                real64 const inputDt,
                                real64 const (&inputGravityVector)[3],
                                string const fluidModelKey );

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geos::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;


    /// The number of jump dofs per element.
    static constexpr int numWdofs = 3;

    /// Constructor.
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            dispEqnRowIndices{ 0 },
      dispColIndices{ 0 },
      jumpEqnRowIndices{ 0 },
      jumpColIndices{ 0 },
      localDispResidual{ 0.0 },
      localJumpResidual{ 0.0 },
      localKww{ { 0.0 } },
      localKwu{ { 0.0 } },
      localKuw{ { 0.0 } },
      localEqMStress { 0.0 },
      localKwpm{ 0.0 },
      localKwpf( 0.0 ),
      wLocal(),
      dispLocal(),
      deltaDispLocal(),
      hInv(),
      xLocal(),
      tractionVec(),
      dTractiondw{ { 0.0 } },
      constitutiveStiffness()
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex jumpEqnRowIndices[numWdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex jumpColIndices[numWdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localDispResidual[numUdofs];

    /// C-array storage for the element local Rw residual vector.
    real64 localJumpResidual[numWdofs];

    /// C-array storage for the element local Kww matrix.
    real64 localKww[numWdofs][numWdofs];

    /// C-array storage for the element local Kwu matrix.
    real64 localKwu[numWdofs][numUdofs];

    /// C-array storage for the element local Kuw matrix.
    real64 localKuw[numUdofs][numWdofs];

    /// C-array storage for the element local EqM*effStress vector.
    real64 localEqMStress[numWdofs];

    /// C-array storage for the element local Kwpm matrix.
    real64 localKwpm[numWdofs];

    /// C-array storage for the element local Kwpf matrix.
    real64 localKwpf;

    /// Stack storage for the element local jump vector
    real64 wLocal[3];

    /// Stack storage for the element displacement vector.
    real64 dispLocal[numUdofs];

    // Stack storage for incremental displacement
    real64 deltaDispLocal[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for Area/Volume
    real64 hInv;

    /// local nodal coordinates
    real64 xLocal[ numNodesPerElem ][ 3 ];

    /// Stack storage for the traction
    real64 tractionVec[3];

    /// Stack storage for the derivative of the traction
    real64 dTractiondw[3][3];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[ 6 ][ 6 ];
  };
  //*****************************************************************************

  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );
  //END_kernelLauncher


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geos::finiteElement::ImplicitKernelBase::setup
   *
   * For the SinglePhasePoromechanicsEFEM implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack,
                              FUNC && kernelOp = NoOpFunc{} ) const;

  /**
   * @copydoc geos::finiteElement::ImplicitKernelBase::complete
   */
  GEOS_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const;

protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_deltaDisp;

  arrayView2d< real64 const > const m_w;

  /// The effective stress at the current time
  arrayView3d< real64 const, solid::STRESS_USD > m_effStress;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_matrixPresDofNumber;

  arrayView1d< globalIndex const > const m_fracturePresDofNumber;

  arrayView1d< globalIndex const > const m_wDofNumber;

  /// The rank global fluid mass
  arrayView1d< real64 const > const m_fluidMass;
  arrayView1d< real64 const > const m_fluidMass_n;
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_dFluidMass;

  /// The rank global densities
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_fluidDensity;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_matrixPressure;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_fracturePressure;

  arrayView2d< real64 const > const m_tractionVec;

  arrayView3d< real64 const > const m_dTraction_dJump;

  arrayView1d< real64 const > const m_dTraction_dPressure;

  arrayView2d< real64 const > const m_nVec;

  arrayView2d< real64 const > const m_tVec1;

  arrayView2d< real64 const > const m_tVec2;

  arrayView2d< real64 const > const m_surfaceCenter;

  arrayView1d< real64 const > const m_surfaceArea;

  arrayView1d< real64 const > const m_elementVolumeCell;

  SortedArrayView< localIndex const > const m_fracturedElems;

  ArrayOfArraysView< localIndex const > const m_cellsToEmbeddedSurfaces;

  /// The gravity vector.
  real64 const m_gravityVector[3];
  real64 const m_gravityAcceleration;

};


using SinglePhaseKernelFactory = finiteElement::KernelFactory< SinglePhasePoromechanicsEFEM,
                                                               EmbeddedSurfaceSubRegion const &,
                                                               arrayView1d< globalIndex const > const,
                                                               arrayView1d< globalIndex const > const,
                                                               string const,
                                                               globalIndex const,
                                                               CRSMatrixView< real64, globalIndex const > const,
                                                               arrayView1d< real64 > const,
                                                               real64 const,
                                                               real64 const (&)[3],
                                                               string const >;

} // namespace poromechanicsEFEMKernels

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEFEM_HPP_
