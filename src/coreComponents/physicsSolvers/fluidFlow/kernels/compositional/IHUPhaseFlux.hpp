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
 * @file IHUPhaseFlux.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHUPHASEFLUX_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHUPHASEFLUX_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernelUtilities
{

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

using Deriv = constitutive::multifluid::DerivativeOffset;

/************************* HELPERS ******************/
namespace UpwindHelpers
{

template< localIndex numComp >
GEOS_HOST_DEVICE
static void assignToZero( real64 & deriv_dP, real64 ( & deriv_dC )[numComp] )
{
  deriv_dP = 0.0;
  for( localIndex ic = 0; ic < numComp; ++ic )
  {
    deriv_dC[ic] = 0.0;
  }
}

template< localIndex numComp >
GEOS_HOST_DEVICE
static void assignToZero( real64 & value, real64 & deriv_dP, real64 ( & deriv_dC )[numComp] )
{
  value = 0.0;
  assignToZero( deriv_dP, deriv_dC );
}

template< localIndex numComp, localIndex numFluxSupportPoints >
GEOS_HOST_DEVICE
static void assignToZero( real64 & value, real64 ( & deriv_dP )[numFluxSupportPoints], real64 ( & deriv_dC )[numFluxSupportPoints][numComp] )
{
  value = 0;
  for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
  {
    assignToZero( deriv_dP[ke], deriv_dC[ke] );
  }
}

template< localIndex numComp >
GEOS_HOST_DEVICE
static void addToDerivativesScaled( real64 const ( &dDeriv_dP ), real64 const ( &dDeriv_dC )[numComp],
                                    real64 const & factor,
                                    real64 ( &deriv_dP ), real64 ( & deriv_dC )[numComp] )
{
  deriv_dP += dDeriv_dP * factor;
  for( localIndex ic = 0; ic < numComp; ++ic )
    deriv_dC[ic] += dDeriv_dC[ic] * factor;
}

template< localIndex numComp, localIndex numFluxSupportPoints >
GEOS_HOST_DEVICE
static void addToDerivativesScaled( real64 const ( &dDeriv_dP )[numFluxSupportPoints], real64 const ( &dDeriv_dC )[numFluxSupportPoints][numComp],
                                    real64 const & factor,
                                    real64 ( & deriv_dP )[numFluxSupportPoints], real64 ( & deriv_dC )[numFluxSupportPoints][numComp] )
{
  for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
  {
    addToDerivativesScaled( dDeriv_dP[ke], dDeriv_dC[ke], factor, deriv_dP[ke], deriv_dC[ke] );
  }
}

template< localIndex numComp >
GEOS_HOST_DEVICE
static void addToDerivatives( real64 const ( &dDeriv_dP ), real64 const ( &dDeriv_dC )[numComp],
                              real64 ( &deriv_dP ), real64 ( & deriv_dC )[numComp] )
{
  addToDerivativesScaled( dDeriv_dP, dDeriv_dC, 1.0, deriv_dP, deriv_dC );
}

template< localIndex numComp, localIndex numFluxSupportPoints >
GEOS_HOST_DEVICE
static void addToDerivatives( real64 const ( &dDeriv_dP )[numFluxSupportPoints], real64 const ( &dDeriv_dC )[numFluxSupportPoints][numComp],
                              real64 ( & deriv_dP )[numFluxSupportPoints], real64 ( & deriv_dC )[numFluxSupportPoints][numComp] )
{
  for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
  {
    addToDerivativesScaled( dDeriv_dP, dDeriv_dC, 1.0, deriv_dP, deriv_dC );
  }
}

template< localIndex numComp >
GEOS_HOST_DEVICE
static void addToValueAndDerivatives( real64 const & dValue, real64 const ( &dDeriv_dP ), real64 const ( &dDeriv_dC )[numComp],
                                      real64 & value, real64 ( &deriv_dP ), real64 ( & deriv_dC )[numComp] )
{
  value += dValue;
  addToDerivatives( dDeriv_dP, dDeriv_dC, deriv_dP, deriv_dC );
}

template< localIndex numComp, localIndex numFluxSupportPoints >
GEOS_HOST_DEVICE
static void addToValueAndDerivatives( real64 const & dValue, real64 const ( &dDeriv_dP )[numFluxSupportPoints], real64 const ( &dDeriv_dC )[numFluxSupportPoints][numComp],
                                      real64 & value, real64 ( & deriv_dP )[numFluxSupportPoints], real64 ( & deriv_dC )[numFluxSupportPoints][numComp] )
{
  value += dValue;
  for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
  {
    addToDerivatives( dDeriv_dP[ke], dDeriv_dC[ke], deriv_dP[ke], deriv_dC[ke] );
  }
}

template< localIndex numComp, localIndex numFluxSupportPoints >
GEOS_HOST_DEVICE
static void assignMobilityAndDerivatives( localIndex const & ip, localIndex const & upwindDir,
                                          localIndex const (&seri)[numFluxSupportPoints],
                                          localIndex const (&sesri)[numFluxSupportPoints],
                                          localIndex const (&sei)[numFluxSupportPoints],
                                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                          real64 & mobility, real64 & dMobility_dP, real64 ( & dMobility_dC )[numComp] )
{
  localIndex const er_up = seri[upwindDir];
  localIndex const esr_up = sesri[upwindDir];
  localIndex const ei_up = sei[upwindDir];

  if( std::fabs( phaseMob[er_up][esr_up][ei_up][ip] ) > 1e-20 )
  {
    mobility += phaseMob[er_up][esr_up][ei_up][ip];
    dMobility_dP += dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
    for( localIndex ic = 0; ic < numComp; ++ic )
    {
      dMobility_dC[ic] += dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dC + ic];
    }
  }
}

template< localIndex numComp, localIndex numFluxSupportPoints >
GEOS_HOST_DEVICE
static void computeFractionalFlowAndDerivatives( localIndex const & k_up, real64 const & mob, real64 const & dMob_dP, real64 const ( &dMob_dC )[numComp],
                                                 real64 const & totMob, real64 const ( &dTotMob_dP )[numFluxSupportPoints], real64 const ( &dTotMob_dC )[numFluxSupportPoints][numComp],
                                                 real64 & fractionalFlow, real64 ( & dFractionalFlow_dP )[numFluxSupportPoints], real64 ( & dFractionalFlow_dC )[numFluxSupportPoints][numComp] )
{
  // guard against no flow region
  // fractional flow too low to let the upstream phase flow
  if( std::fabs( mob ) > 1e-20 )
  {
    real64 const invTotMob = 1 / totMob;

    fractionalFlow = mob * invTotMob;

    addToDerivativesScaled( dMob_dP, dMob_dC, invTotMob, dFractionalFlow_dP[k_up], dFractionalFlow_dC[k_up] );

    addToDerivativesScaled( dTotMob_dP, dTotMob_dC, -fractionalFlow * invTotMob, dFractionalFlow_dP, dFractionalFlow_dC );
  }
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
upwindMobilityViscous( localIndex const numPhase,
                       localIndex const ip,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const (&transmissibility)[2],
                       real64 const (&dTrans_dPres)[2],
                       real64 const totFlux, // in fine should be a ElemnetViewConst once seq form are in place
                       ElementViewConst< arrayView1d< real64 const > > const & pres,
                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                       ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                       ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                       integer const hasCapPressure,
                       localIndex & upwindDir,
                       real64 & mobility,
                       real64 (&dMobility_dP),
                       real64 (& dMobility_dC)[numComp] )
{
  UPWIND scheme;
  scheme.template getUpwindDirectionViscous< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                      ip,
                                                                                      seri,
                                                                                      sesri,
                                                                                      sei,
                                                                                      transmissibility,
                                                                                      dTrans_dPres,
                                                                                      totFlux,
                                                                                      pres,
                                                                                      gravCoef,
                                                                                      phaseMob,
                                                                                      dCompFrac_dCompDens,
                                                                                      phaseMassDens,
                                                                                      dPhaseMassDens,
                                                                                      dPhaseVolFrac,
                                                                                      phaseCapPressure,
                                                                                      dPhaseCapPressure_dPhaseVolFrac,
                                                                                      hasCapPressure,
                                                                                      upwindDir );

  //reinit
  assignToZero( mobility, dMobility_dP, dMobility_dC );
  assignMobilityAndDerivatives( ip, upwindDir, seri, sesri, sei, phaseMob, dPhaseMob, mobility, dMobility_dP, dMobility_dC );
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
upwindMobilityGravity( localIndex const numPhase,
                       localIndex const ip,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const (&transmissibility)[2],
                       real64 const (&dTrans_dPres)[2],
                       real64 const totFlux,            //in fine should be a ElemnetViewConst once seq form are in place
                       ElementViewConst< arrayView1d< real64 const > > const & pres,
                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                       ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                       ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                       integer const hasCapPressure,
                       integer const checkPhasePresenceInGravity,
                       localIndex & upwindDir,
                       real64 & mobility,
                       real64 ( &dMobility_dP ),
                       real64 ( & dMobility_dC )[numComp] )
{
  UPWIND scheme;
  scheme.template getUpwindDirectionGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                      ip,
                                                                                      seri,
                                                                                      sesri,
                                                                                      sei,
                                                                                      transmissibility,
                                                                                      dTrans_dPres,
                                                                                      totFlux,
                                                                                      pres,
                                                                                      gravCoef,
                                                                                      phaseMob,
                                                                                      dCompFrac_dCompDens,
                                                                                      phaseMassDens,
                                                                                      dPhaseMassDens,
                                                                                      phaseVolFrac,
                                                                                      dPhaseVolFrac,
                                                                                      phaseCapPressure,
                                                                                      dPhaseCapPressure_dPhaseVolFrac,
                                                                                      hasCapPressure,
                                                                                      checkPhasePresenceInGravity,
                                                                                      upwindDir );

  //reinit
  assignToZero( mobility, dMobility_dP, dMobility_dC );
  assignMobilityAndDerivatives( ip, upwindDir, seri, sesri, sei, phaseMob, dPhaseMob, mobility, dMobility_dP, dMobility_dC );
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
upwindMobilityCapillary( localIndex const numPhase,
                         localIndex const ip,
                         localIndex const (&seri)[numFluxSupportPoints],
                         localIndex const (&sesri)[numFluxSupportPoints],
                         localIndex const (&sei)[numFluxSupportPoints],
                         real64 const (&transmissibility)[2],
                         real64 const (&dTrans_dPres)[2],
                         real64 const totFlux,          //in fine should be a ElemnetViewConst once seq form are in place
                         ElementViewConst< arrayView1d< real64 const > > const & pres,
                         ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                         ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                         ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                         ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                         ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                         integer const hasCapPressure,
                         localIndex & upwindDir,
                         real64 & mobility,
                         real64 ( &dMobility_dP ),
                         real64 ( & dMobility_dC )[numComp] )
{
  UPWIND scheme;
  scheme.template getUpwindDirectionCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                        ip,
                                                                                        seri,
                                                                                        sesri,
                                                                                        sei,
                                                                                        transmissibility,
                                                                                        dTrans_dPres,
                                                                                        totFlux,
                                                                                        pres,
                                                                                        gravCoef,
                                                                                        phaseMob,
                                                                                        dCompFrac_dCompDens,
                                                                                        phaseMassDens,
                                                                                        dPhaseMassDens,
                                                                                        dPhaseVolFrac,
                                                                                        phaseCapPressure,
                                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                                        hasCapPressure,
                                                                                        upwindDir );

  //reinit
  assignToZero( mobility, dMobility_dP, dMobility_dC );
  assignMobilityAndDerivatives( ip, upwindDir, seri, sesri, sei, phaseMob, dPhaseMob, mobility, dMobility_dP, dMobility_dC );
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
computeFractionalFlowViscous( localIndex const numPhase,
                              localIndex const ip,
                              localIndex const (&seri)[numFluxSupportPoints],
                              localIndex const (&sesri)[numFluxSupportPoints],
                              localIndex const (&sei)[numFluxSupportPoints],
                              real64 const (&transmissibility)[2],
                              real64 const (&dTrans_dPres)[2],
                              real64 const totFlux,
                              real64 const totMob,
                              real64 const (&dTotMob_dP)[numFluxSupportPoints],
                              real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
                              ElementViewConst< arrayView1d< real64 const > > const & pres,
                              ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                              ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                              ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                              ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                              ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                              ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                              integer const hasCapPressure,
                              real64 & fractionalFlow,
                              real64 (& dFractionalFlow_dP)[numFluxSupportPoints],
                              real64 (& dFractionalFlow_dC)[numFluxSupportPoints][numComp] )
{
  localIndex k_up;
  real64 mob{};
  real64 dMob_dP{};
  real64 dMob_dC[numComp]{};

  upwindMobilityViscous< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                  ip,
                                                                  seri,
                                                                  sesri,
                                                                  sei,
                                                                  transmissibility,
                                                                  dTrans_dPres,
                                                                  totFlux,
                                                                  pres,
                                                                  gravCoef,
                                                                  dCompFrac_dCompDens,
                                                                  phaseMassDens,
                                                                  dPhaseMassDens,
                                                                  phaseMob,
                                                                  dPhaseMob,
                                                                  dPhaseVolFrac,
                                                                  phaseCapPressure,
                                                                  dPhaseCapPressure_dPhaseVolFrac,
                                                                  hasCapPressure,
                                                                  k_up,
                                                                  mob,
                                                                  dMob_dP,
                                                                  dMob_dC );

  // reinit
  assignToZero( fractionalFlow, dFractionalFlow_dP, dFractionalFlow_dC );
  computeFractionalFlowAndDerivatives( k_up, mob, dMob_dP, dMob_dC,
                                       totMob, dTotMob_dP, dTotMob_dC,
                                       fractionalFlow, dFractionalFlow_dP, dFractionalFlow_dC );
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
computeFractionalFlowGravity( localIndex const numPhase,
                              localIndex const ip,
                              localIndex const (&seri)[numFluxSupportPoints],
                              localIndex const (&sesri)[numFluxSupportPoints],
                              localIndex const (&sei)[numFluxSupportPoints],
                              real64 const (&transmissibility)[2],
                              real64 const (&dTrans_dPres)[2],
                              real64 const totFlux,
                              real64 const totMob,
                              real64 const (&dTotMob_dP)[numFluxSupportPoints],
                              real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
                              ElementViewConst< arrayView1d< real64 const > > const & pres,
                              ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                              ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                              ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                              ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                              ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                              ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                              ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                              integer const hasCapPressure,
                              integer const checkPhasePresenceInGravity,
                              real64 & fractionalFlow,
                              real64 (& dFractionalFlow_dP)[numFluxSupportPoints],
                              real64 (& dFractionalFlow_dC)[numFluxSupportPoints][numComp] )
{
  localIndex k_up;
  real64 mob{};
  real64 dMob_dP{};
  real64 dMob_dC[numComp]{};

  upwindMobilityGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                  ip,
                                                                  seri,
                                                                  sesri,
                                                                  sei,
                                                                  transmissibility,
                                                                  dTrans_dPres,
                                                                  totFlux,
                                                                  pres,
                                                                  gravCoef,
                                                                  dCompFrac_dCompDens,
                                                                  phaseMassDens,
                                                                  dPhaseMassDens,
                                                                  phaseMob,
                                                                  dPhaseMob,
                                                                  phaseVolFrac,
                                                                  dPhaseVolFrac,
                                                                  phaseCapPressure,
                                                                  dPhaseCapPressure_dPhaseVolFrac,
                                                                  hasCapPressure,
                                                                  checkPhasePresenceInGravity,
                                                                  k_up,
                                                                  mob,
                                                                  dMob_dP,
                                                                  dMob_dC );

  // reinit
  assignToZero( fractionalFlow, dFractionalFlow_dP, dFractionalFlow_dC );
  computeFractionalFlowAndDerivatives( k_up, mob, dMob_dP, dMob_dC,
                                       totMob, dTotMob_dP, dTotMob_dC,
                                       fractionalFlow, dFractionalFlow_dP, dFractionalFlow_dC );
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
computeFractionalFlowCapillary( localIndex const numPhase,
                                localIndex const ip,
                                localIndex const (&seri)[numFluxSupportPoints],
                                localIndex const (&sesri)[numFluxSupportPoints],
                                localIndex const (&sei)[numFluxSupportPoints],
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const totFlux,
                                real64 const totMob,
                                real64 const (&dTotMob_dP)[numFluxSupportPoints],
                                real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
                                ElementViewConst< arrayView1d< real64 const > > const & pres,
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                integer const hasCapPressure,
                                real64 & fractionalFlow,
                                real64 ( & dFractionalFlow_dP)[numFluxSupportPoints],
                                real64 ( & dFractionalFlow_dC)[numFluxSupportPoints][numComp]
                                )
{
  localIndex k_up;
  real64 mob{};
  real64 dMob_dP{};
  real64 dMob_dC[numComp]{};

  upwindMobilityCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                    ip,
                                                                    seri,
                                                                    sesri,
                                                                    sei,
                                                                    transmissibility,
                                                                    dTrans_dPres,
                                                                    totFlux,
                                                                    pres,
                                                                    gravCoef,
                                                                    dCompFrac_dCompDens,
                                                                    phaseMassDens,
                                                                    dPhaseMassDens,
                                                                    phaseMob,
                                                                    dPhaseMob,
                                                                    dPhaseVolFrac,
                                                                    phaseCapPressure,
                                                                    dPhaseCapPressure_dPhaseVolFrac,
                                                                    hasCapPressure,
                                                                    k_up,
                                                                    mob,
                                                                    dMob_dP,
                                                                    dMob_dC );

  // reinit
  assignToZero( fractionalFlow, dFractionalFlow_dP, dFractionalFlow_dC );
  computeFractionalFlowAndDerivatives( k_up, mob, dMob_dP, dMob_dC,
                                       totMob, dTotMob_dP, dTotMob_dC,
                                       fractionalFlow, dFractionalFlow_dP, dFractionalFlow_dC );
}

/**
 * @brief  Struct defining formation of potential from different Physics (flagged by enum type T) to be used
 *            in Upwind discretization schemes
 * @tparam numComp number of component
 * @tparam T the concerned physics (Viscou,Gravity or Capillary)
 * @tparam numFluxSupportPoints number of point in the stencil
 */
struct computePotentialViscous
{
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void compute( localIndex const GEOS_UNUSED_PARAM( numPhase ),
                       localIndex const GEOS_UNUSED_PARAM( ip ),
                       localIndex const (&GEOS_UNUSED_PARAM( seri ))[numFluxSupportPoints],
                       localIndex const (&GEOS_UNUSED_PARAM( sesri ))[numFluxSupportPoints],
                       localIndex const (&GEOS_UNUSED_PARAM( sei ))[numFluxSupportPoints],
                       real64 const (&GEOS_UNUSED_PARAM( transmissibility ))[2],
                       real64 const (&GEOS_UNUSED_PARAM( dTrans_dPres ))[2],
                       real64 const totFlux,
                       ElementViewConst< arrayView1d< real64 const > > const & GEOS_UNUSED_PARAM( gravCoef ),
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const &
                       GEOS_UNUSED_PARAM( dCompFrac_dCompDens ),
                       ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const &
                       GEOS_UNUSED_PARAM( phaseMassDens ),
                       ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const &
                       GEOS_UNUSED_PARAM( dPhaseMassDens ),
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const &
                       GEOS_UNUSED_PARAM( dPhaseVolFrac ),
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const &
                       GEOS_UNUSED_PARAM( phaseCapPressure ),
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const &
                       GEOS_UNUSED_PARAM( dPhaseCapPressure_dPhaseVolFrac ),
                       real64 & pot,
                       real64( &GEOS_UNUSED_PARAM( dPot_dPres ))[numFluxSupportPoints],
                       real64( &GEOS_UNUSED_PARAM( dPot_dComp ))[numFluxSupportPoints][numComp] )
  {
    pot = totFlux;
    //could be relevant for symmetry to include derivative
  }
};

/*! @copydoc computePotential
 */
struct computePotentialGravity
{
  /*! @copydoc computePotential::compute
   *
   * @brief specialization for gravitational driving forces which only relies on total flux
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void compute( localIndex const GEOS_UNUSED_PARAM( numPhase ),
                       localIndex const ip,
                       integer const checkPhasePresenceInGravity,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const (&transmissibility)[2],
                       real64 const (&dTrans_dPres)[2],
                       real64 const GEOS_UNUSED_PARAM( totFlux ),
                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                       ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                       ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const &
                       GEOS_UNUSED_PARAM( dPhaseVolFrac ),
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const &
                       GEOS_UNUSED_PARAM( phaseCapPressure ),
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const &
                       GEOS_UNUSED_PARAM( dPhaseCapPressure_dPhaseVolFrac ),
                       real64 & pot,
                       real64 ( & dPot_dPres )[numFluxSupportPoints],
                       real64 ( & dPot_dComp )[numFluxSupportPoints][numComp] )
  {
    //init
    assignToZero( pot, dPot_dPres, dPot_dComp );

    //working arrays
    real64 densMean{};
    real64 dDensMean_dPres[numFluxSupportPoints]{};
    real64 dDensMean_dComp[numFluxSupportPoints][numComp]{};
    isothermalCompositionalMultiphaseFVMKernels::helpers::
      calculateMeanDensity( ip, seri, sesri, sei,
                            checkPhasePresenceInGravity,
                            phaseVolFrac, dCompFrac_dCompDens,
                            phaseMassDens, dPhaseMassDens,
                            densMean, dDensMean_dPres, dDensMean_dComp );

    // compute potential difference MPFA-style
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      real64 const gravD = transmissibility[i] * gravCoef[er][esr][ei];
      real64 const dGravD_dP = dTrans_dPres[i] * gravCoef[er][esr][ei];
      pot += densMean * gravD;
      dPot_dPres[i] += densMean * dGravD_dP;

      // need to add contributions from both cells the mean density depends on
      addToDerivativesScaled( dDensMean_dPres, dDensMean_dComp, gravD, dPot_dPres, dPot_dComp );
    }
  }

};

/*! @copydoc computePotential
 */
struct computePotentialCapillary
{
  /*! @copydoc computePotential::compute
   *
   * @brief specialization for capillary driving forces which only relies on total flux
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void compute( localIndex const numPhase,
                       localIndex const ip,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const (&transmissibility)[2],
                       real64 const (&dTrans_dPres)[2],
                       real64 const GEOS_UNUSED_PARAM( totFlux ),
                       ElementViewConst< arrayView1d< real64 const > > const & GEOS_UNUSED_PARAM( gravCoef ),
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const &
                       GEOS_UNUSED_PARAM( dCompFrac_dCompDens ),
                       ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const &
                       GEOS_UNUSED_PARAM( phaseMassDens ),
                       ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const &
                       GEOS_UNUSED_PARAM( dPhaseMassDens ),
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                       real64 & pot,
                       real64 ( & dPot_dPres)[numFluxSupportPoints],
                       real64 (& dPot_dComp)[numFluxSupportPoints][numComp] )
  {
    //init
    assignToZero( pot, dPot_dPres, dPot_dComp );

    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      real64 const capPres = phaseCapPressure[er][esr][ei][0][ip];

      pot += transmissibility[i] * capPres;

      dPot_dPres[i] += dTrans_dPres[i] * capPres;

      // need to add contributions from both cells
      for( localIndex jp = 0; jp < numPhase; ++jp )
      {
        real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];

        dPot_dPres[i] += transmissibility[i] * dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dP];

        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPot_dComp[i][jc] += transmissibility[i] * dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dC + jc];
        }
      }
    }
  }
};

/// Form potential-related parts of fluxes

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void computePotentialFluxesGravity( localIndex const numPhase,
                                           localIndex const ip,
                                           localIndex const (&seri)[numFluxSupportPoints],
                                           localIndex const (&sesri)[numFluxSupportPoints],
                                           localIndex const (&sei)[numFluxSupportPoints],
                                           real64 const (&transmissibility)[2],
                                           real64 const (&dTrans_dPres)[2],
                                           real64 const totFlux,
                                           real64 const totMob,
                                           real64 const (&dTotMob_dP)[numFluxSupportPoints],
                                           real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
                                           ElementViewConst< arrayView1d< real64 const > > const & pres,
                                           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                                           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                           ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                           ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                           localIndex const hasCapPressure,
                                           integer const checkPhasePresenceInGravity,
                                           real64 & phaseFlux,
                                           real64 (& dPhaseFlux_dP)[numFluxSupportPoints],
                                           real64 ( & dPhaseFlux_dC)[numFluxSupportPoints][numComp] )
{

  real64 fflow{};
  real64 dFflow_dP[numFluxSupportPoints]{};
  real64 dFflow_dC[numFluxSupportPoints][numComp]{};

  real64 pot{};
  real64 dPot_dP[numFluxSupportPoints]{};
  real64 dPot_dC[numFluxSupportPoints][numComp]{};

  //
  UpwindHelpers::computePotentialGravity::compute< numComp, numFluxSupportPoints >( numPhase,
                                                                                    ip,
                                                                                    checkPhasePresenceInGravity,
                                                                                    seri,
                                                                                    sesri,
                                                                                    sei,
                                                                                    transmissibility,
                                                                                    dTrans_dPres,
                                                                                    totFlux,
                                                                                    gravCoef,
                                                                                    dCompFrac_dCompDens,
                                                                                    phaseMassDens,
                                                                                    dPhaseMassDens,
                                                                                    phaseVolFrac,
                                                                                    dPhaseVolFrac,
                                                                                    phaseCapPressure,
                                                                                    dPhaseCapPressure_dPhaseVolFrac,
                                                                                    pot,
                                                                                    dPot_dP,
                                                                                    dPot_dC );

  // and the fractional flow for gravitational part as \lambda_i^{up}/\sum_{numPhase}(\lambda_k^{up}) with up decided upon
  // the Upwind strategy
  UpwindHelpers::computeFractionalFlowGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                        ip,
                                                                                        seri,
                                                                                        sesri,
                                                                                        sei,
                                                                                        transmissibility,
                                                                                        dTrans_dPres,
                                                                                        totFlux,
                                                                                        totMob,
                                                                                        dTotMob_dP,
                                                                                        dTotMob_dC,
                                                                                        pres,
                                                                                        gravCoef,
                                                                                        dCompFrac_dCompDens,
                                                                                        phaseMassDens,
                                                                                        dPhaseMassDens,
                                                                                        phaseMob,
                                                                                        dPhaseMob,
                                                                                        phaseVolFrac,
                                                                                        dPhaseVolFrac,
                                                                                        phaseCapPressure,
                                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                                        hasCapPressure,
                                                                                        checkPhasePresenceInGravity,
                                                                                        fflow,
                                                                                        dFflow_dP,
                                                                                        dFflow_dC );


  for( localIndex jp = 0; jp < numPhase; ++jp )
  {
    if( ip != jp )
    {

      real64 potOther{};
      real64 dPotOther_dP[numFluxSupportPoints]{};
      real64 dPotOther_dC[numFluxSupportPoints][numComp]{};

      //Fetch pot for phase j!=i defined as \rho_j g dz/dx
      UpwindHelpers::computePotentialGravity::compute< numComp, numFluxSupportPoints >( numPhase,
                                                                                        jp,
                                                                                        checkPhasePresenceInGravity,
                                                                                        seri,
                                                                                        sesri,
                                                                                        sei,
                                                                                        transmissibility,
                                                                                        dTrans_dPres,
                                                                                        totFlux,
                                                                                        gravCoef,
                                                                                        dCompFrac_dCompDens,
                                                                                        phaseMassDens,
                                                                                        dPhaseMassDens,
                                                                                        phaseVolFrac,
                                                                                        dPhaseVolFrac,
                                                                                        phaseCapPressure,
                                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                                        potOther,
                                                                                        dPotOther_dP,
                                                                                        dPotOther_dC );

      //Eventually get the mobility of the second phase
      localIndex k_up_o = -1;
      real64 mobOther{};
      real64 dMobOther_dP{};
      real64 dMobOther_dC[numComp]{};

      // and the other mobility for gravitational part as \lambda_j^{up} with up decided upon
      // the Upwind strategy - Note that it should be the same as the gravitational fractional flow

      UpwindHelpers::upwindMobilityGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                     jp,
                                                                                     seri,
                                                                                     sesri,
                                                                                     sei,
                                                                                     transmissibility,
                                                                                     dTrans_dPres,
                                                                                     totFlux,
                                                                                     pres,
                                                                                     gravCoef,
                                                                                     dCompFrac_dCompDens,
                                                                                     phaseMassDens,
                                                                                     dPhaseMassDens,
                                                                                     phaseMob,
                                                                                     dPhaseMob,
                                                                                     phaseVolFrac,
                                                                                     dPhaseVolFrac,
                                                                                     phaseCapPressure,
                                                                                     dPhaseCapPressure_dPhaseVolFrac,
                                                                                     hasCapPressure,
                                                                                     checkPhasePresenceInGravity,
                                                                                     k_up_o,
                                                                                     mobOther,
                                                                                     dMobOther_dP,
                                                                                     dMobOther_dC );


      // Assembling gravitational flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (G_i - G_k)
      phaseFlux -= fflow * mobOther * (pot - potOther);

      // mobOther part
      UpwindHelpers::addToDerivativesScaled( dMobOther_dP, dMobOther_dC, -fflow * (pot - potOther), dPhaseFlux_dP[k_up_o], dPhaseFlux_dC[k_up_o] );

      // mob related part of dFflow_dP is only upstream defined but totMob related is defined everywhere
      UpwindHelpers::addToDerivativesScaled( dFflow_dP, dFflow_dC, -mobOther * (pot - potOther), dPhaseFlux_dP, dPhaseFlux_dC );

      // potential difference part
      UpwindHelpers::addToDerivativesScaled( dPot_dP, dPot_dC, -fflow * mobOther, dPhaseFlux_dP, dPhaseFlux_dC );
      UpwindHelpers::addToDerivativesScaled( dPotOther_dP, dPotOther_dC, fflow * mobOther, dPhaseFlux_dP, dPhaseFlux_dC );
    }
  }

}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void computePotentialFluxesCapillary( localIndex const numPhase,
                                             localIndex const ip,
                                             localIndex const (&seri)[numFluxSupportPoints],
                                             localIndex const (&sesri)[numFluxSupportPoints],
                                             localIndex const (&sei)[numFluxSupportPoints],
                                             real64 const (&transmissibility)[2],
                                             real64 const (&dTrans_dPres)[2],
                                             real64 const totFlux,
                                             real64 const totMob,
                                             real64 const (&dTotMob_dP)[numFluxSupportPoints],
                                             real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
                                             ElementViewConst< arrayView1d< real64 const > > const & pres,
                                             ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                             ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                             ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                             ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                             ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                             ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                             ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                             ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                             ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                             localIndex const hasCapPressure,
                                             real64 & phaseFlux,
                                             real64 (& dPhaseFlux_dP)[numFluxSupportPoints],
                                             real64 ( & dPhaseFlux_dC)[numFluxSupportPoints][numComp] )
{
  real64 fflow{};
  real64 dFflow_dP[numFluxSupportPoints]{};
  real64 dFflow_dC[numFluxSupportPoints][numComp]{};

  real64 pot{};
  real64 dPot_dP[numFluxSupportPoints]{};
  real64 dPot_dC[numFluxSupportPoints][numComp]{};

  computePotentialCapillary::compute< numComp, numFluxSupportPoints >( numPhase,
                                                                       ip,
                                                                       seri,
                                                                       sesri,
                                                                       sei,
                                                                       transmissibility,
                                                                       dTrans_dPres,
                                                                       totFlux,
                                                                       gravCoef,
                                                                       dCompFrac_dCompDens,
                                                                       phaseMassDens,
                                                                       dPhaseMassDens,
                                                                       dPhaseVolFrac,
                                                                       phaseCapPressure,
                                                                       dPhaseCapPressure_dPhaseVolFrac,
                                                                       pot,
                                                                       dPot_dP,
                                                                       dPot_dC );

  // and the fractional flow for gravitational part as \lambda_i^{up}/\sum_{numPhase}(\lambda_k^{up}) with up decided upon
  // the Upwind strategy
  UpwindHelpers::computeFractionalFlowCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                          ip,
                                                                                          seri,
                                                                                          sesri,
                                                                                          sei,
                                                                                          transmissibility,
                                                                                          dTrans_dPres,
                                                                                          totFlux,
                                                                                          totMob,
                                                                                          dTotMob_dP,
                                                                                          dTotMob_dC,
                                                                                          pres,
                                                                                          gravCoef,
                                                                                          dCompFrac_dCompDens,
                                                                                          phaseMassDens,
                                                                                          dPhaseMassDens,
                                                                                          phaseMob,
                                                                                          dPhaseMob,
                                                                                          dPhaseVolFrac,
                                                                                          phaseCapPressure,
                                                                                          dPhaseCapPressure_dPhaseVolFrac,
                                                                                          hasCapPressure,
                                                                                          fflow,
                                                                                          dFflow_dP,
                                                                                          dFflow_dC );


  for( localIndex jp = 0; jp < numPhase; ++jp )
  {
    if( ip != jp )
    {
      real64 potOther{};
      real64 dPotOther_dP[numFluxSupportPoints]{};
      real64 dPotOther_dC[numFluxSupportPoints][numComp]{};

      //Fetch pot for phase j!=i defined as \rho_j g dz/dx
      computePotentialCapillary::compute< numComp, numFluxSupportPoints >( numPhase,
                                                                           jp,
                                                                           seri,
                                                                           sesri,
                                                                           sei,
                                                                           transmissibility,
                                                                           dTrans_dPres,
                                                                           totFlux,
                                                                           gravCoef,
                                                                           dCompFrac_dCompDens,
                                                                           phaseMassDens,
                                                                           dPhaseMassDens,
                                                                           dPhaseVolFrac,
                                                                           phaseCapPressure,
                                                                           dPhaseCapPressure_dPhaseVolFrac,
                                                                           potOther,
                                                                           dPotOther_dP,
                                                                           dPotOther_dC );

      //Eventually get the mobility of the second phase
      localIndex k_up_o = -1;
      real64 mobOther{};
      real64 dMobOther_dP{};
      real64 dMobOther_dC[numComp]{};

      // and the other mobility for gravitational part as \lambda_j^{up} with up decided upon
      // the Upwind strategy - Note that it should be the same as the gravitational fractional flow

      UpwindHelpers::upwindMobilityCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                       jp,
                                                                                       seri,
                                                                                       sesri,
                                                                                       sei,
                                                                                       transmissibility,
                                                                                       dTrans_dPres,
                                                                                       totFlux,
                                                                                       pres,
                                                                                       gravCoef,
                                                                                       dCompFrac_dCompDens,
                                                                                       phaseMassDens,
                                                                                       dPhaseMassDens,
                                                                                       phaseMob,
                                                                                       dPhaseMob,
                                                                                       dPhaseVolFrac,
                                                                                       phaseCapPressure,
                                                                                       dPhaseCapPressure_dPhaseVolFrac,
                                                                                       hasCapPressure,
                                                                                       k_up_o,
                                                                                       mobOther,
                                                                                       dMobOther_dP,
                                                                                       dMobOther_dC );


      // Assembling gravitational flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (G_i - G_k)
      phaseFlux -= fflow * mobOther * (pot - potOther);

      // mobOther part
      UpwindHelpers::addToDerivativesScaled( dMobOther_dP, dMobOther_dC, -fflow * (pot - potOther), dPhaseFlux_dP[k_up_o], dPhaseFlux_dC[k_up_o] );

      // mob related part of dFflow_dP is only upstream defined but totMob related is defined everywhere
      UpwindHelpers::addToDerivativesScaled( dFflow_dP, dFflow_dC, -mobOther * (pot - potOther), dPhaseFlux_dP, dPhaseFlux_dC );

      // potential difference part
      UpwindHelpers::addToDerivativesScaled( dPot_dP, dPot_dC, -fflow * mobOther, dPhaseFlux_dP, dPhaseFlux_dC );
      UpwindHelpers::addToDerivativesScaled( dPotOther_dP, dPotOther_dC, fflow * mobOther, dPhaseFlux_dP, dPhaseFlux_dC );
    }
  }
}

}//end of struct UpwindHelpers

/************************* UPWIND ******************/

/**
 * @brief Template base class for different upwind Scheme
 * @tparam T physics concerned by the scheme if specialized
 */
class UpwindScheme
{

public:

  //default ctor
  UpwindScheme() = default;

  //usual copy ctor
  UpwindScheme( UpwindScheme const & scheme ) = default;

  //default move ctor
  UpwindScheme( UpwindScheme && ) = default;

  //deleted copy and move assignement
  UpwindScheme & operator=( UpwindScheme const & ) = delete;

  UpwindScheme & operator=( UpwindScheme && ) = delete;

  virtual ~UpwindScheme() = default;

  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOS_HOST_DEVICE
  inline
  void getUpwindDirectionViscous( localIndex const numPhase,
                                  localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  real64 const (&transmissibility)[2],
                                  real64 const (&dTrans_dPres)[2],
                                  real64 const totFlux,            //in fine should be a ElemnetViewConst once seq form are in place
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                  integer const hasCapPressure,
                                  localIndex & upwindDir
                                  )
  {
    real64 pot{};
    /// each derived concrete class has to define a computePotential method that is calling UpwindScheme::potential method with a specific
    /// lamda defining how to get these potentials
    UPWIND::template computePotentialViscous< numComp, numFluxSupportPoints >( numPhase,
                                                                               ip,
                                                                               seri,
                                                                               sesri,
                                                                               sei,
                                                                               transmissibility,
                                                                               dTrans_dPres,
                                                                               totFlux,
                                                                               pres,
                                                                               gravCoef,
                                                                               phaseMob,
                                                                               dCompFrac_dCompDens,
                                                                               phaseMassDens,
                                                                               dPhaseMassDens,
                                                                               dPhaseVolFrac,
                                                                               phaseCapPressure,
                                                                               dPhaseCapPressure_dPhaseVolFrac,
                                                                               hasCapPressure,
                                                                               pot );
    //all definition has been changed to fit pot>0 => first cell is upstream
    upwindDir = (pot > 0) ? 0 : 1;
  }


  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOS_HOST_DEVICE
  void getUpwindDirectionGravity( localIndex const numPhase,
                                  localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  real64 const (&transmissibility)[2],
                                  real64 const (&dTrans_dPres)[2],
                                  real64 const totFlux,
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                  integer const hasCapPressure,
                                  integer const checkPhasePresenceInGravity,
                                  localIndex & upwindDir
                                  )
  {
    real64 pot{};

    /// each derived concrete class has to define a computePotential method that is calling UpwindScheme::potential method with a specific
    /// lamda defining how to get these potentials
    UPWIND::template computePotentialGravity< numComp, numFluxSupportPoints >( numPhase,
                                                                               ip,
                                                                               seri,
                                                                               sesri,
                                                                               sei,
                                                                               transmissibility,
                                                                               dTrans_dPres,
                                                                               totFlux,
                                                                               pres,
                                                                               gravCoef,
                                                                               phaseMob,
                                                                               dCompFrac_dCompDens,
                                                                               phaseMassDens,
                                                                               dPhaseMassDens,
                                                                               phaseVolFrac,
                                                                               dPhaseVolFrac,
                                                                               phaseCapPressure,
                                                                               dPhaseCapPressure_dPhaseVolFrac,
                                                                               hasCapPressure,
                                                                               checkPhasePresenceInGravity,
                                                                               pot );

    //all definition has been changed to fit pot>0 => first cell is upstream
    upwindDir = (pot >= 0) ? 0 : 1;
  }


  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOS_HOST_DEVICE
  void getUpwindDirectionCapillary( localIndex const numPhase,
                                    localIndex const ip,
                                    localIndex const (&seri)[numFluxSupportPoints],
                                    localIndex const (&sesri)[numFluxSupportPoints],
                                    localIndex const (&sei)[numFluxSupportPoints],
                                    real64 const (&transmissibility)[2],
                                    real64 const (&dTrans_dPres)[2],
                                    real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
                                    ElementViewConst< arrayView1d< real64 const > > const & pres,
                                    ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                    ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                    ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                    ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                    ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                    ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                    ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                    ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                    integer const hasCapPressure,
                                    localIndex & upwindDir
                                    )
  {
    real64 pot{};

    // each derived concrete class has to define a computePotential method that is calling UpwindScheme::potential method with a specific
    // lamda defining how to get these potentials
    UPWIND::template computePotentialCapillary< numComp, numFluxSupportPoints >( numPhase,
                                                                                 ip,
                                                                                 seri,
                                                                                 sesri,
                                                                                 sei,
                                                                                 transmissibility,
                                                                                 dTrans_dPres,
                                                                                 totFlux,
                                                                                 pres,
                                                                                 gravCoef,
                                                                                 phaseMob,
                                                                                 dCompFrac_dCompDens,
                                                                                 phaseMassDens,
                                                                                 dPhaseMassDens,
                                                                                 dPhaseVolFrac,
                                                                                 phaseCapPressure,
                                                                                 dPhaseCapPressure_dPhaseVolFrac,
                                                                                 hasCapPressure,
                                                                                 pot );

    //all definition has been changed to fit pot>0 => first cell is upstream
    upwindDir = (pot >= 0) ? 0 : 1;
  }



  // templated way of evaluating the potential (to the exception of viscous one) which relies on
  // up-or-downwinded mobility terms pre-multiplying potential differences
  template< localIndex numComp, localIndex numFluxSupportPoints, typename LAMBDA >
  GEOS_HOST_DEVICE
  static void potential( localIndex numPhase,
                         localIndex ip,
                         localIndex const (&seri)[numFluxSupportPoints],
                         localIndex const (&sesri)[numFluxSupportPoints],
                         localIndex const (&sei)[numFluxSupportPoints],
                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                         real64 & weightedPotential,
                         LAMBDA && fn )
  {
    //getPhase Pot
    real64 pot{};
    real64 pot_dP[numFluxSupportPoints]{};
    real64 pot_dC[numFluxSupportPoints][numComp]{};

    fn( ip, pot, pot_dP, pot_dC );

    localIndex const k_up = 0;
    localIndex const k_dw = 1;

    //loop other other phases to form
    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( jp != ip )
      {
        localIndex const er_up = seri[k_up];
        localIndex const esr_up = sesri[k_up];
        localIndex const ei_up = sei[k_up];

        localIndex const er_dw = seri[k_dw];
        localIndex const esr_dw = sesri[k_dw];
        localIndex const ei_dw = sei[k_dw];

        real64 potOther{};
        real64 potOther_dP[numFluxSupportPoints]{};
        real64 potOther_dC[numFluxSupportPoints][numComp]{};

        fn( jp, potOther, potOther_dP, potOther_dC );

        real64 const mob_up = phaseMob[er_up][esr_up][ei_up][jp];
        real64 const mob_dw = phaseMob[er_dw][esr_dw][ei_dw][jp];

        weightedPotential += (pot - potOther >= 0) ? mob_dw * (potOther - pot) : mob_up * (potOther - pot);

      }
    }
  }

};

/**
 * @brief  Class describing the Hybrid Upwind scheme as defined in "Consistent upwinding for sequential fully implicit
 *         multiscale compositional simulation" (Moncorge,2020)
 */
class HybridUpwind : public UpwindScheme
{

public:
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static
  void computePotentialViscous( localIndex const numPhase,
                                localIndex const ip,
                                localIndex const (&seri)[numFluxSupportPoints],
                                localIndex const (&sesri)[numFluxSupportPoints],
                                localIndex const (&sei)[numFluxSupportPoints],
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const totalFlux,
                                ElementViewConst< arrayView1d< real64 const > > const & GEOS_UNUSED_PARAM( pres ),
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const &
                                GEOS_UNUSED_PARAM( phaseMob ),
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                integer const GEOS_UNUSED_PARAM( hasCapPressure ),
                                real64 & potential )
  {
    real64 dPot_dP[numFluxSupportPoints]{};
    real64 dPot_dC[numFluxSupportPoints][numComp]{};
    UpwindHelpers::computePotentialViscous::compute< numComp, numFluxSupportPoints >(
      numPhase,
      ip,
      seri,
      sesri,
      sei,
      transmissibility,
      dTrans_dPres,
      totalFlux,
      gravCoef,
      dCompFrac_dCompDens,
      phaseMassDens,
      dPhaseMassDens,
      dPhaseVolFrac,
      phaseCapPressure,
      dPhaseCapPressure_dPhaseVolFrac,
      potential,
      dPot_dP,
      dPot_dC );
  }

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static
  void computePotentialGravity( localIndex const numPhase,
                                localIndex const ip,
                                localIndex const (&seri)[numFluxSupportPoints],
                                localIndex const (&sesri)[numFluxSupportPoints],
                                localIndex const (&sei)[numFluxSupportPoints],
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const totalFlux,
                                ElementViewConst< arrayView1d< real64 const > > const & GEOS_UNUSED_PARAM( pres ),
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                integer const GEOS_UNUSED_PARAM( hasCapPressure ),
                                integer const checkPhasePresenceInGravity,
                                real64 & potential )
  {

    //Form total velocity
    potential = 0;

    //the arg lambda allows us to access some genericity
    UpwindScheme::template potential< numComp, numFluxSupportPoints >( numPhase, ip, seri, sesri, sei,
                                                                       phaseMob, potential,
                                                                       [&]( localIndex ipp,
                                                                            real64 & potential_,
                                                                            real64 (& dPotential_dP_)[numFluxSupportPoints],
                                                                            real64 (& dPotential_dC_)[numFluxSupportPoints][numComp] ) {

      UpwindHelpers::computePotentialGravity::compute< numComp, numFluxSupportPoints >(
        numPhase,
        ipp,
        checkPhasePresenceInGravity,
        seri,
        sesri,
        sei,
        transmissibility,
        dTrans_dPres,
        totalFlux,
        gravCoef,
        dCompFrac_dCompDens,
        phaseMassDens,
        dPhaseMassDens,
        phaseVolFrac,
        dPhaseVolFrac,
        phaseCapPressure,
        dPhaseCapPressure_dPhaseVolFrac,
        potential_,
        dPotential_dP_,
        dPotential_dC_ );

    } );
  }


  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static
  void computePotentialCapillary( localIndex const numPhase,
                                  localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  real64 const (&transmissibility)[2],
                                  real64 const (&dTrans_dPres)[2],
                                  real64 const totalFlux,
                                  ElementViewConst< arrayView1d< real64 const > > const & GEOS_UNUSED_PARAM( pres ),
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const &
                                  phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                  integer const GEOS_UNUSED_PARAM( hasCapPressure ),
                                  real64 & potential )

  {

    //Form total velocity
    potential = 0;

    //the arg lambda allows us to access some genericity
    UpwindScheme::template potential< numComp, numFluxSupportPoints >( numPhase, ip, seri, sesri, sei,
                                                                       phaseMob, potential,
                                                                       [&]( localIndex ipp,
                                                                            real64 & potential_,
                                                                            real64 (& dPotential_dP_)[numFluxSupportPoints],
                                                                            real64 (& dPotential_dC_)[numFluxSupportPoints][numComp] )
    {

      UpwindHelpers::computePotentialCapillary::compute< numComp, numFluxSupportPoints >(
        numPhase,
        ipp,
        seri,
        sesri,
        sei,
        transmissibility,
        dTrans_dPres,
        totalFlux,
        gravCoef,
        dCompFrac_dCompDens,
        phaseMassDens,
        dPhaseMassDens,
        dPhaseVolFrac,
        phaseCapPressure,
        dPhaseCapPressure_dPhaseVolFrac,
        potential_,
        dPotential_dP_,
        dPotential_dC_ );
    } );
  }

};

/*** IHU ***/

struct IHUPhaseFlux
{

  using UPWIND_SCHEME = HybridUpwind;

  static constexpr double minTotMob = 1e-12;

  /**
   * @brief Form the Implicit Hybrid Upwind from pressure gradient and gravitational head
   * @tparam numComp number of components
   * @tparam numFluxSupportPoints number of flux support points
   * @param numPhase number of phases
   * @param ip phase index
   * @param hasCapPressure flag indicating if there is capillary pressure
   * @param seri arraySlice of the stencil-implied element region index
   * @param sesri arraySlice of the stencil-implied element subregion index
   * @param sei arraySlice of the stencil-implied element index
   * @param trans transmissibility at the connection
   * @param dTrans_dPres derivative of transmissibility wrt pressure
   * @param pres pressure
   * @param gravCoef gravitational coefficient
   * @param phaseMob phase mobility
   * @param dPhaseMob derivative of phase mobility wrt pressure, temperature, comp density
   * @param dPhaseVolFrac derivative of phase volume fraction wrt pressure, temperature, comp density
   * @param dCompFrac_dCompDens derivative of component fraction wrt component density
   * @param phaseMassDens phase mass density
   * @param dPhaseMassDens derivative of phase mass density wrt pressure, temperature, comp fraction
   * @param phaseCapPressure phase capillary pressure
   * @param dPhaseCapPressure_dPhaseVolFrac derivative of phase capillary pressure wrt phase volume fraction
   * @param potGrad potential gradient for this phase
   * @param phaseFlux phase flux
   * @param dPhaseFlux_dP derivative of phase flux wrt pressure
   * @param dPhaseFlux_dC derivative of phase flux wrt comp density
   */
  template< integer numComp, integer numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute( integer const numPhase,
           integer const ip,
           integer const hasCapPressure,
           integer const checkPhasePresenceInGravity,
           localIndex const ( &seri )[numFluxSupportPoints],
           localIndex const ( &sesri )[numFluxSupportPoints],
           localIndex const ( &sei )[numFluxSupportPoints],
           real64 const ( &trans )[2],
           real64 const ( &dTrans_dPres )[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
           ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
           ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
           real64 & potGrad,
           real64 ( &phaseFlux ),
           real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
           real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp] )
  {

    //loop over all phases to form total velocity
    real64 totFlux{};
    real64 dTotFlux_dP[numFluxSupportPoints]{};
    real64 dTotFlux_dC[numFluxSupportPoints][numComp]{};

    //store totMob upwinded by PPU for later schemes
    real64 totMob{};
    real64 dTotMob_dP[numFluxSupportPoints]{};
    real64 dTotMob_dC[numFluxSupportPoints][numComp]{};

    //unelegant but need dummy when forming PPU total velocity
    real64 dPhaseFlux_dTrans;

    for( integer jp = 0; jp < numPhase; ++jp )
    {
      PPUPhaseFlux::compute( numPhase, jp,
                             hasCapPressure, checkPhasePresenceInGravity,
                             seri, sesri, sei,
                             trans, dTrans_dPres,
                             pres, gravCoef,
                             phaseMob, dPhaseMob,
                             phaseVolFrac, dPhaseVolFrac,
                             dCompFrac_dCompDens,
                             phaseMassDens, dPhaseMassDens,
                             phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac,
                             potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, dPhaseFlux_dTrans );

      // accumulate into total flux
      UpwindHelpers::addToValueAndDerivatives( phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC,
                                               totFlux, dTotFlux_dP, dTotFlux_dC );
/*
      // old wrong way:
      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        totMob += phaseMob[seri[ke]][sesri[ke]][sei[ke]][jp];
        dTotMob_dP[ke] += dPhaseMob[seri[ke]][sesri[ke]][sei[ke]][jp][Deriv::dP];
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dTotMob_dC[ke][jc] += dPhaseMob[seri[ke]][sesri[ke]][sei[ke]][jp][Deriv::dC + jc];
        }
      }
 */

      // new correct way:
      // choose upstream cell
      localIndex const k_up = (phaseFlux >= 0) ? 0 : 1;
      // accumulate into total mobility
      UpwindHelpers::assignMobilityAndDerivatives( jp, k_up, seri, sesri, sei,
                                                   phaseMob, dPhaseMob,
                                                   totMob, dTotMob_dP[k_up], dTotMob_dC[k_up] );
    }

    // safeguard
    totMob = LvArray::math::max( totMob, minTotMob );

    // assign to zero
    UpwindHelpers::assignToZero( phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

    // fractional flow loop with IHU
    // maybe needed to have density out for upwinding

    // choose upstream cell
    // create local work arrays
    real64 viscousPhaseFlux{};
    real64 dViscousPhaseFlux_dP[numFluxSupportPoints]{};
    real64 dViscousPhaseFlux_dC[numFluxSupportPoints][numComp]{};

    real64 fractionalFlow{};
    real64 dFractionalFlow_dP[numFluxSupportPoints]{};
    real64 dFractionalFlow_dC[numFluxSupportPoints][numComp]{};

    // and the fractional flow for viscous part as \lambda_i^{up}/\sum_{NP}(\lambda_j^{up}) with up decided upon
    // the Upwind strategy
    UpwindHelpers::computeFractionalFlowViscous< numComp, numFluxSupportPoints,
                                                 UPWIND_SCHEME >( numPhase,
                                                                  ip,
                                                                  seri,
                                                                  sesri,
                                                                  sei,
                                                                  trans,
                                                                  dTrans_dPres,
                                                                  totFlux,
                                                                  totMob,
                                                                  dTotMob_dP,
                                                                  dTotMob_dC,
                                                                  pres,
                                                                  gravCoef,
                                                                  dCompFrac_dCompDens,
                                                                  phaseMassDens,
                                                                  dPhaseMassDens,
                                                                  phaseMob,
                                                                  dPhaseMob,
                                                                  dPhaseVolFrac,
                                                                  phaseCapPressure,
                                                                  dPhaseCapPressure_dPhaseVolFrac,
                                                                  hasCapPressure,
                                                                  fractionalFlow,
                                                                  dFractionalFlow_dP,
                                                                  dFractionalFlow_dC );

    /// Assembling the viscous flux (and derivatives) from fractional flow and total velocity as \phi_{\mu} = f_i^{up,\mu} uT
    viscousPhaseFlux = fractionalFlow * totFlux;

    UpwindHelpers::addToDerivativesScaled( dFractionalFlow_dP, dFractionalFlow_dC,
                                           totFlux,
                                           dViscousPhaseFlux_dP, dViscousPhaseFlux_dC );

    //NON-FIXED UT -- to be canceled out if considered fixed
    UpwindHelpers::addToDerivativesScaled( dTotFlux_dP, dTotFlux_dC,
                                           fractionalFlow,
                                           dViscousPhaseFlux_dP, dViscousPhaseFlux_dC );

    // accumulate in the flux and its derivatives
    UpwindHelpers::addToValueAndDerivatives( viscousPhaseFlux, dViscousPhaseFlux_dP, dViscousPhaseFlux_dC,
                                             phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

    /// Assembling the gravitational flux (and derivatives) from fractional flow and total velocity as \phi_{g} = f_i^{up,g} uT

    real64 gravitationalPhaseFlux{};
    real64 gravitationalPhaseFlux_dP[numFluxSupportPoints]{};
    real64 gravitationalPhaseFlux_dC[numFluxSupportPoints][numComp]{};

    UpwindHelpers::computePotentialFluxesGravity< numComp,
                                                  numFluxSupportPoints, UPWIND_SCHEME >(
      numPhase,
      ip,
      seri,
      sesri,
      sei,
      trans,
      dTrans_dPres,
      totFlux,
      totMob,
      dTotMob_dP,
      dTotMob_dC,
      pres,
      gravCoef,
      phaseMob,
      dPhaseMob,
      phaseVolFrac,
      dPhaseVolFrac,
      dCompFrac_dCompDens,
      phaseMassDens,
      dPhaseMassDens,
      phaseCapPressure,
      dPhaseCapPressure_dPhaseVolFrac,
      hasCapPressure,
      checkPhasePresenceInGravity,
      gravitationalPhaseFlux,
      gravitationalPhaseFlux_dP,
      gravitationalPhaseFlux_dC );

    //update phaseFlux from gravitational
    UpwindHelpers::addToValueAndDerivatives( gravitationalPhaseFlux, gravitationalPhaseFlux_dP, gravitationalPhaseFlux_dC,
                                             phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

    if( hasCapPressure )
    {
      /// Assembling the capillary flux (and derivatives) from fractional flow and total velocity as \phi_{g} = f_i^{up,g} uT

      real64 capillaryPhaseFlux{};
      real64 capillaryPhaseFlux_dP[numFluxSupportPoints]{};
      real64 capillaryPhaseFlux_dC[numFluxSupportPoints][numComp]{};

      UpwindHelpers::computePotentialFluxesCapillary< numComp,
                                                      numFluxSupportPoints, UPWIND_SCHEME >(
        numPhase,
        ip,
        seri,
        sesri,
        sei,
        trans,
        dTrans_dPres,
        totFlux,
        totMob,
        dTotMob_dP,
        dTotMob_dC,
        pres,
        gravCoef,
        phaseMob,
        dPhaseMob,
        dPhaseVolFrac,
        dCompFrac_dCompDens,
        phaseMassDens,
        dPhaseMassDens,
        phaseCapPressure,
        dPhaseCapPressure_dPhaseVolFrac,
        hasCapPressure,
        capillaryPhaseFlux,
        capillaryPhaseFlux_dP,
        capillaryPhaseFlux_dC );

      //update phaseFlux from capillary
      UpwindHelpers::addToValueAndDerivatives( capillaryPhaseFlux, capillaryPhaseFlux_dP, capillaryPhaseFlux_dC,
                                               phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

    }//end if cappres
  }

};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geos


#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHUPHASEFLUX_HPP
