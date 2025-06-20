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
 * @file CriticalVolume.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_CRITICALVOLUME_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_CRITICALVOLUME_HPP_

#include "ModelParameters.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class CriticalVolume : public ModelParameters
{
public:
  CriticalVolume( std::unique_ptr< ModelParameters > parameters );
  ~CriticalVolume() override = default;

  static std::unique_ptr< ModelParameters > create( std::unique_ptr< ModelParameters > parameters );

  array1d< real64 > m_componentCriticalVolume;

  struct viewKeyStruct
  {
    static constexpr char const * componentCriticalVolumeString() { return "componentCriticalVolume"; }
  };

protected:
  void registerParametersImpl( MultiFluidBase * fluid ) override;

  void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;

  /**
   * @brief Estimate critical volumes using Ihmels' (2010) correlation
   * @details reference: http://dx.doi.org/10.1021/je100167w
   * @param[in] numComponents The number of components
   * @param[in] criticalPressure The component critical pressures
   * @param[in] criticalTemperature The component critical temperatures
   * @param[in] criticalVolume The component critical volumes
   */
  static void calculateCriticalVolume( integer const numComponents,
                                       arrayView1d< const real64 > const criticalPressure,
                                       arrayView1d< const real64 > const criticalTemperature,
                                       arrayView1d< real64 > const criticalVolume );
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_CRITICALVOLUME_HPP_
