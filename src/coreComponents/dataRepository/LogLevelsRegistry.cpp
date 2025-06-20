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

#include "LogLevelsRegistry.hpp"
#include <unordered_set>

namespace geos
{

void LogLevelsRegistry::addEntry( integer condition, std::string_view description )
{

  auto & targetValues = m_logLevelsDescriptions[condition];

  if( !(std::find( targetValues.begin(), targetValues.end(), description ) != targetValues.end()))
  {
    targetValues.emplace_back( std::string( description ) );
  }
}

string LogLevelsRegistry::buildLogLevelDescription() const
{
  std::ostringstream description;
  description << "Sets the level of information to write in the standard output (the console typically).\n"
                 "Information output from lower logLevels is added with the desired log level";
  for( auto const & [logLevel, logDescriptions] : m_logLevelsDescriptions )
  {
    description << GEOS_FMT( "\n{}\n", logLevel );
    for( size_t idxDesc = 0; idxDesc < logDescriptions.size(); idxDesc++ )
    {
      description << " - " << logDescriptions[idxDesc];
      if( idxDesc != logDescriptions.size() - 1 )
        description << '\n';
    }
  }
  return description.str();
}

}
