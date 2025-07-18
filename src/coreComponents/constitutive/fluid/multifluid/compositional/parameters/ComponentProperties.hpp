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
 * @file ComponentProperties.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_COMPONENTPROPERTIES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_COMPONENTPROPERTIES_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/**
 * @brief Class holding standard component properties for a compositional fluid model.
 */
class ComponentProperties final
{
public:
  ComponentProperties( string_array const & componentNames,
                       array1d< real64 > const & componentMolarWeight );
  ~ComponentProperties() = default;
  ComponentProperties( const ComponentProperties & ) = default;
  const ComponentProperties & operator=( const ComponentProperties & ) = delete;

  /**
   * @brief Get the number of components
   * @return The number of components
   */
  integer getNumberOfComponents() const { return m_componentNames.size(); }

  /**
   * @brief Assign component types to the components
   */
  void classifyComponents()
  {
    classifyComponents( m_componentNames, m_componentType );
  }

  /**
   * Data accessors
   */
  string_array const & getComponentName() const { return m_componentNames; }
  arrayView1d< integer > const & getComponentType() const { return m_componentType; }
  arrayView1d< real64 > const & getComponentMolarWeight() const { return m_componentMolarWeight; }
  arrayView1d< real64 > const & getComponentCriticalPressure() const { return m_componentCriticalPressure; }
  arrayView1d< real64 > const & getComponentCriticalTemperature() const { return m_componentCriticalTemperature; }
  arrayView1d< real64 > const & getComponentAcentricFactor() const { return m_componentAcentricFactor; }
  arrayView1d< real64 > const & getComponentVolumeShift() const { return m_componentVolumeShift; }

  struct KernelWrapper
  {
    KernelWrapper( arrayView1d< integer const > const & componentType,
                   arrayView1d< real64 const > const & componentMolarWeight,
                   arrayView1d< real64 const > const & componentCriticalPressure,
                   arrayView1d< real64 const > const & componentCriticalTemperature,
                   arrayView1d< real64 const > const & componentAcentricFactor,
                   arrayView1d< real64 const > const & componentVolumeShift,
                   arrayView2d< real64 const > const & componentBinaryCoeff );

    /**
     * @brief Move the KernelWrapper to the given execution space, optionally touching it.
     * @param space the space to move the KernelWrapper to
     * @param touch whether the KernelWrapper should be touched in the new space or not
     * @note This function exists to enable holding KernelWrapper objects in an ArrayView
     *       and have their contents properly moved between memory spaces.
     */
    void move( LvArray::MemorySpace const space, bool const touch )
    {
      m_componentType.move( space, touch );
      m_componentMolarWeight.move( space, touch );
      m_componentCriticalPressure.move( space, touch );
      m_componentCriticalTemperature.move( space, touch );
      m_componentAcentricFactor.move( space, touch );
      m_componentVolumeShift.move( space, touch );
      m_componentBinaryCoeff.move( space, touch );
    }

    // Standard compositional input
    arrayView1d< integer const > m_componentType;
    arrayView1d< real64 const > m_componentMolarWeight;
    arrayView1d< real64 const > m_componentCriticalPressure;
    arrayView1d< real64 const > m_componentCriticalTemperature;
    arrayView1d< real64 const > m_componentAcentricFactor;
    arrayView1d< real64 const > m_componentVolumeShift;
    arrayView2d< real64 const > m_componentBinaryCoeff;
  };

  /**
   * @brief Function to create and return a KernelWrapper
   * @return the KernelWrapper object
   */
  KernelWrapper createKernelWrapper() const;

private:
  static void classifyComponents( string_array const & componentNames, array1d< integer > & componentType );

public:
  // Standard compositional input
  string_array const & m_componentNames;
  array1d< real64 > const & m_componentMolarWeight;
  array1d< integer > m_componentType;
  array1d< real64 > m_componentCriticalPressure;
  array1d< real64 > m_componentCriticalTemperature;
  array1d< real64 > m_componentAcentricFactor;
  array1d< real64 > m_componentVolumeShift;
  array2d< real64 > m_componentBinaryCoeff;
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_COMPONENTPROPERTIES_HPP_
