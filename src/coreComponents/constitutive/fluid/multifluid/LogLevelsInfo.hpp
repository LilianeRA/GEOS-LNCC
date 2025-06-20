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
 * @file LogLevelsInfo.hpp
 * This file contains common log level informations for fluid
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_LOGLEVELSINFO_HPP
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_LOGLEVELSINFO_HPP

#include "common/DataTypes.hpp"

namespace geos
{

namespace logInfo
{
/**
 * @name Common LogLevels info structures. They must comply with the `is_log_level_info` trait.
 */
///@{

struct PVT
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Output PVT table"; }
};

/// @cond DO_NOT_DOCUMENT
/// @endcond
///@}

}

}

#endif // GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_LOGLEVELSINFO_HPP
