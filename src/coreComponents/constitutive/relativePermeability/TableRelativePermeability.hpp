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
 * @file TableRelativePermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP
#define GEOS_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityInterpolators.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{
namespace constitutive
{

class TableRelativePermeability : public RelativePermeabilityBase
{
public:

  /// order of the phase properties for three-phase flow
  struct ThreePhasePairPhaseType
  {
    enum : integer
    {
      WETTING = 0,                ///< wetting phase property
      INTERMEDIATE_WETTING = 1,   ///< intermediate phase property
      NONWETTING = 2,             ///< non-wetting phase property
      INTERMEDIATE_NONWETTING = 3 ///< intermediate phase property
    };
  };

  /// order of the phase properties in the wetting-non-wetting data
  struct TwoPhasePairPhaseType
  {
    enum : integer
    {
      WETTING = 0,   ///< wetting phase property
      NONWETTING = 1 ///< non-wetting phase property
    };
  };

  TableRelativePermeability( std::string const & name, dataRepository::Group * const parent );

  static std::string catalogName() { return "TableRelativePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }


  /// Type of kernel wrapper for in-kernel update
  class KernelWrapper final : public RelativePermeabilityBaseUpdate
  {
public:

    KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & relPermKernelWrappers,
                   arrayView1d< real64 const > const & phaseMinVolumeFraction,
                   real64 const & waterOilPhaseMaxVolumeFraction,
                   arrayView1d< integer const > const & phaseTypes,
                   arrayView1d< integer const > const & phaseOrder,
                   ThreePhaseInterpolator const & threePhaseInterpolator,
                   arrayView3d< real64, constitutive::relperm::USD_RELPERM > const & phaseRelPerm,
                   arrayView4d< real64, constitutive::relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
                   arrayView3d< real64, constitutive::relperm::USD_RELPERM > const & phaseTrappedVolFrac );


    GEOS_HOST_DEVICE
    void computeTwoPhase( integer const ipWetting,
                          integer const ipNonWetting,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                          arraySlice1d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                          arraySlice2d< real64, constitutive::relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOS_HOST_DEVICE
    void computeThreePhase( integer const ipWetting,
                            integer const ipInter,
                            integer const ipNonWetting,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                            arraySlice1d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                            arraySlice2d< real64, constitutive::relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOS_HOST_DEVICE
    void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                  arraySlice1d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseTrappedVolFrac,
                  arraySlice1d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                  arraySlice2d< real64, constitutive::relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOS_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override;

private:

    /// Kernel wrappers for relative permeabilities in the following order:
    /// Two-phase flow:
    ///  0- wetting-phase
    ///  1- non-wetting-phase
    /// Three-phase flow:
    ///  0- wetting-phase
    ///  1- intermediate phase (wetting-intermediate data)
    ///  2- non-wetting-phase
    ///  3- intermediate phase (non-wetting-intermediate data)
    arrayView1d< TableFunction::KernelWrapper const > m_relPermKernelWrappers;

    /// Minimum volume fraction for each phase (deduced from the table)
    arrayView1d< real64 const > m_phaseMinVolumeFraction;

    real64 const m_waterOilRelPermMaxValue;

    ThreePhaseInterpolator const m_threePhaseInterpolator;
  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr char const * relPermKernelWrappersString() { return "relPermWrappers"; }
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
    static constexpr char const * waterOilMaxRelPermString() { return "waterOilMaxRelPerm"; }
    static constexpr char const * wettingNonWettingRelPermTableNamesString() { return "wettingNonWettingRelPermTableNames"; }
    static constexpr char const * wettingIntermediateRelPermTableNamesString() { return "wettingIntermediateRelPermTableNames"; }
    static constexpr char const * nonWettingIntermediateRelPermTableNamesString() { return "nonWettingIntermediateRelPermTableNames"; }
    static constexpr char const * threePhaseInterpolatorString() { return "threePhaseInterpolator"; }
  };

  arrayView1d< real64 const > getPhaseMinVolumeFraction() const override { return m_phaseMinVolumeFraction; };

private:

  virtual void postInputInitialization() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  /// Relative permeability table names (one for each phase in the wetting-non-wetting pair)
  string_array m_wettingNonWettingRelPermTableNames;

  /// Relative permeability table names (one for each phase in the wetting-intermediate pair)
  string_array m_wettingIntermediateRelPermTableNames;

  /// Relative permeability table names (one for each phase in the non-wetting-intermediate pair)
  string_array m_nonWettingIntermediateRelPermTableNames;

  /// Kernel wrappers for relative permeabilities in the following order:
  /// Two-phase flow:
  ///  0- wetting-phase
  ///  1- non-wetting-phase
  /// Three-phase flow:
  ///  0- wetting-phase
  ///  1- intermediate phase (wetting-intermediate data)
  ///  2- non-wetting-phase
  ///  3- intermediate phase (non-wetting-intermediate data)
  array1d< TableFunction::KernelWrapper > m_relPermKernelWrappers;

  /// Min phase volume fractions (deduced from the tables). With Baker, only the water phase entry is used
  array1d< real64 > m_phaseMinVolumeFraction;

  real64 m_waterOilMaxRelPerm;

  ThreePhaseInterpolator m_threePhaseInterpolator;

};

GEOS_HOST_DEVICE
inline void
TableRelativePermeability::KernelWrapper::
  computeTwoPhase( integer const ipWetting,
                   integer const ipNonWetting,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                   arraySlice1d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                   arraySlice2d< real64, constitutive::relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  using TPT = TableRelativePermeability::TwoPhasePairPhaseType;

  // water rel perm
  phaseRelPerm[ipWetting] =
    m_relPermKernelWrappers[TPT::WETTING].compute( &(phaseVolFraction)[ipWetting],
                                                   &(dPhaseRelPerm_dPhaseVolFrac)[ipWetting][ipWetting] );

  // oil rel perm
  phaseRelPerm[ipNonWetting] =
    m_relPermKernelWrappers[TPT::NONWETTING].compute( &(phaseVolFraction)[ipNonWetting],
                                                      &(dPhaseRelPerm_dPhaseVolFrac)[ipNonWetting][ipNonWetting] );

}



GEOS_HOST_DEVICE
inline void
TableRelativePermeability::KernelWrapper::
  computeThreePhase( integer const ipWetting,
                     integer const ipInter,
                     integer const ipNonWetting,
                     arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                     arraySlice1d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                     arraySlice2d< real64, constitutive::relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  real64 interRelPerm_wi = 0; // oil rel perm using two-phase gas-oil data
  real64 dInterRelPerm_wi_dInterVolFrac = 0; // derivative w.r.t to So
  real64 interRelPerm_nwi = 0; // oil rel perm using two-phase gas-oil data
  real64 dInterRelPerm_nwi_dInterVolFrac = 0; // derivative w.r.t to So
  using TPT = TableRelativePermeability::ThreePhasePairPhaseType;

  // 1) Wetting and intermediate phase relative permeabilities using two-phase wetting-intermediate data

  // wetting rel perm
  phaseRelPerm[ipWetting] =
    m_relPermKernelWrappers[TPT::WETTING].compute( &(phaseVolFraction)[ipWetting],
                                                   &(dPhaseRelPerm_dPhaseVolFrac)[ipWetting][ipWetting] );

  // intermediate rel perm
  interRelPerm_wi =
    m_relPermKernelWrappers[TPT::INTERMEDIATE_WETTING].compute( &(phaseVolFraction)[ipInter],
                                                                &dInterRelPerm_wi_dInterVolFrac );


  // 2) Non-wetting and intermediate phase relative permeabilities using two-phase non-wetting-intermediate data

  // gas rel perm
  phaseRelPerm[ipNonWetting] =
    m_relPermKernelWrappers[TPT::NONWETTING].compute( &(phaseVolFraction)[ipNonWetting],
                                                      &(dPhaseRelPerm_dPhaseVolFrac)[ipNonWetting][ipNonWetting] );

  // oil rel perm
  interRelPerm_nwi =
    m_relPermKernelWrappers[TPT::INTERMEDIATE_NONWETTING].compute( &(phaseVolFraction)[ipInter],
                                                                   &dInterRelPerm_nwi_dInterVolFrac );

  // 3) Compute the "three-phase" oil relperm

  // use saturation-weighted interpolation
  real64 const shiftedWettingVolFrac = (phaseVolFraction[ipWetting] - m_phaseMinVolumeFraction[ipWetting]);

  if( m_threePhaseInterpolator == ThreePhaseInterpolator::BAKER )
  {
    relpermInterpolators::Baker::compute( shiftedWettingVolFrac,
                                          phaseVolFraction[ipNonWetting],
                                          m_phaseOrder,
                                          interRelPerm_wi,
                                          dInterRelPerm_wi_dInterVolFrac,
                                          interRelPerm_nwi,
                                          dInterRelPerm_nwi_dInterVolFrac,
                                          phaseRelPerm[ipInter],
                                          dPhaseRelPerm_dPhaseVolFrac[ipInter] );

  }
  else// if( m_threePhaseInterpolator == ThreePhaseInterpolator::STONEII )
  {
    relpermInterpolators::Stone2::compute( shiftedWettingVolFrac,
                                           phaseVolFraction[ipNonWetting],
                                           m_phaseOrder,
                                           m_waterOilRelPermMaxValue,
                                           interRelPerm_wi,
                                           dInterRelPerm_wi_dInterVolFrac,
                                           interRelPerm_nwi,
                                           dInterRelPerm_nwi_dInterVolFrac,
                                           phaseRelPerm[ipWetting],
                                           dPhaseRelPerm_dPhaseVolFrac[ipWetting][ipWetting],
                                           phaseRelPerm[ipNonWetting],
                                           dPhaseRelPerm_dPhaseVolFrac[ipNonWetting][ipNonWetting],
                                           phaseRelPerm[ipInter],
                                           dPhaseRelPerm_dPhaseVolFrac[ipInter] );
  }



}


GEOS_HOST_DEVICE
inline void
TableRelativePermeability::KernelWrapper::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
           arraySlice1d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseTrappedVolFrac,
           arraySlice1d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseRelPerm,
           arraySlice2d< real64, constitutive::relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseRelPerm_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  using PT = RelativePermeabilityBase::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];
  integer const ipOil   = m_phaseOrder[PT::OIL];
  integer const ipGas   = m_phaseOrder[PT::GAS];

  if( ipWater >= 0 && ipOil >= 0 && ipGas >= 0 )
  {
    computeThreePhase( ipWater, // wetting
                       ipOil, // intermediate
                       ipGas, // non-wetting
                       phaseVolFraction,
                       phaseRelPerm,
                       dPhaseRelPerm_dPhaseVolFrac );

  }
  else if( ipWater < 0 )
  {
    computeTwoPhase( ipOil, // wetting
                     ipGas, // non-wetting
                     phaseVolFraction,
                     phaseRelPerm,
                     dPhaseRelPerm_dPhaseVolFrac );
  }
  else if( ipOil < 0 )
  {
    computeTwoPhase( ipWater, // wetting
                     ipGas, // non-wetting
                     phaseVolFraction,
                     phaseRelPerm,
                     dPhaseRelPerm_dPhaseVolFrac );
  }
  else if( ipGas < 0 )
  {
    computeTwoPhase( ipWater, // wetting
                     ipOil, // non-wetting
                     phaseVolFraction,
                     phaseRelPerm,
                     dPhaseRelPerm_dPhaseVolFrac );
  }

  // update trapped phase volume fraction
  if( ipWater >= 0 )
  {
    phaseTrappedVolFrac[ipWater] = LvArray::math::min( phaseVolFraction[ipWater], m_phaseMinVolumeFraction[ipWater] );
  }
  if( ipGas >= 0 )
  {
    phaseTrappedVolFrac[ipGas] = LvArray::math::min( phaseVolFraction[ipGas], m_phaseMinVolumeFraction[ipGas] );
  }
  if( ipOil >= 0 )
  {
    phaseTrappedVolFrac[ipOil] = LvArray::math::min( phaseVolFraction[ipOil], m_phaseMinVolumeFraction[ipOil] );
  }

}

GEOS_HOST_DEVICE
inline void
TableRelativePermeability::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          arraySlice1d< geos::real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const
{
  compute( phaseVolFraction,
           m_phaseTrappedVolFrac[k][q],
           m_phaseRelPerm[k][q],
           m_dPhaseRelPerm_dPhaseVolFrac[k][q] );

}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP
