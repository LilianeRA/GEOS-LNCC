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
 * @file StringUtilities.hpp
 */

#ifndef GEOS_COMMON_FORMAT_STRINGUTILITIES_HPP_
#define GEOS_COMMON_FORMAT_STRINGUTILITIES_HPP_

#include "common/DataTypes.hpp"

#include <iomanip>
#include <sstream>

namespace geos
{
namespace stringutilities
{

/**
 * @brief Return a copy of the string in lower case.
 * @param input The input string which is not modified.
 * @return A new string instance.
 */
string toLower( string const & input );

/**
 * @brief Join strings or other printable objects with a delimiter.
 * @tparam IT    type of iterator into the range of objects to join
 * @tparam S     type of delimiter, usually char, char const * or string
 * @param first  iterator to start of the range
 * @param last   iterator past-the-end of the range
 * @param delim  delimiter used to glue together strings
 * @return a string containing input values concatenated with a delimiter
 */
template< typename IT, typename S = char >
string join( IT first, IT last, S const & delim = S() )
{
  if( first == last )
  {
    return {};
  }
  std::ostringstream oss;
  oss << *first;
  while( ++first != last )
  {
    oss << delim << *first;
  }
  return oss.str();
}

/**
 * @brief Join strings or other printable objects with a delimiter.
 * @tparam CONTAINER type of container to join
 * @tparam S the type of delimiter, usually char, char const * or string
 * @param container the container to join
 * @param delim delimiter used to glue together strings
 * @return a string containing input values concatenated with a delimiter
 */
template< typename CONTAINER, typename S = char >
string join( CONTAINER const & container, S const & delim = S() )
{
  return join( std::begin( container ), std::end( container ), delim );
}

/**
 * @brief Join strings or other printable objects returned by a formatter functor.
 * @tparam IT                type of iterator into the range of objects to join
 * @tparam S                 type of delimiter, usually char, char const * or string
 * @tparam LAMBDA            type of formatter functor, usually `[]( auto it ) -> string`
 * @param formattingFunc     formatter function to get a formattable value from a IT iterator
 * @param first              iterator to start of the range
 * @param last               iterator past-the-end of the range
 * @param delim              delimiter used to glue together strings
 * @return a string containing input values concatenated with a delimiter
 */
template< typename IT, typename S, typename LAMBDA >
string joinLambda( IT first, IT last, S const & delim, LAMBDA formattingFunc )
{
  if( first == last )
  {
    return {};
  }
  std::ostringstream oss;
  oss << formattingFunc( first );
  while( ++first != last )
  {
    oss << delim << formattingFunc( first );
  }
  return oss.str();
}

/**
 * @brief Join strings or other printable objects returned by a formatter functor.
 * @tparam CONTAINER         type of container to join
 * @tparam S                 type of delimiter, usually char, char const * or string
 * @tparam LAMBDA            type of formatter functor, usually `[]( auto it ) -> string`
 * @param formattingFunc     formatter function to get a formattable value from an iterator of the container
 * @param container          container to join
 * @param delim              delimiter used to glue together strings
 * @return a string containing input values concatenated with a delimiter
 */
template< typename CONTAINER, typename S, typename LAMBDA >
string joinLambda( CONTAINER const & container, S const & delim, LAMBDA formattingFunc )
{
  return joinLambda( std::begin( container ), std::end( container ), delim, formattingFunc );
}

/**
 * @brief Concatenate variadic arguments into a string with a delimiter.
 * @tparam S type of delimiter (printable to std::ostringstream)
 * @tparam T type of first argument (printable to std::ostringstream)
 * @tparam Ts types of remaining arguments (printable to std::ostringstream)
 * @param delim delimiter
 * @param v first value
 * @param vs remaining values
 * @return string containing concatenated printed arguments
 */
template< typename S = char, typename T, typename ... Ts >
string concat( S const & delim, T const & v, Ts const & ... vs )
{
  std::ostringstream oss;
  oss << v;
  // Use array initializer and comma trick to get "fold expression" pre C++-17
  using expander = int[];
  (void) expander{ 0, ( void ( oss << delim << vs ), 0) ... };
  return oss.str();
}

/**
 * @brief Subdivide the string in substrings by the specified delimiters.
 * @tparam CONTAINER The templated class of the results container (stdVector by default).
 * @param str The string to subdivide.
 * @param delimiters String that contains the list of possible delimiters.
 * @param treatConsecutiveDelimAsOne If enabled, consecutive delimiters will be treated as one.
 *                                   If not enabled, consecutive delimiters will result in empty entries.
 * @param preTrimStr If enabled, delimiters at the borders of the string will be ignored.
 *                   If not enabled, those delimiters will result in in empty entries.
 * @return The container of the divided substrings.
 */
template< template< class ... > class CONTAINER = stdVector >
CONTAINER< string > tokenize( string const & str,
                              string const & delimiters,
                              bool const treatConsecutiveDelimAsOne = true,
                              bool const preTrimStr = false )
{
  CONTAINER< string > tokens;
  string::size_type tokenBegin, tokenEnd, strEnd;

  if( preTrimStr )
  {
    tokenBegin = str.find_first_not_of( delimiters );
    strEnd = str.find_last_not_of( delimiters ) + 1;
  }
  else
  {
    tokenBegin = 0;
    strEnd = str.size();
  }

  while( ( ( tokenEnd = str.find_first_of( delimiters, tokenBegin ) ) < strEnd ) && tokenBegin < strEnd )
  {
    tokens.emplace_back( str.substr( tokenBegin, tokenEnd - tokenBegin ) );
    tokenBegin = !treatConsecutiveDelimAsOne ? tokenEnd + 1 : str.find_first_not_of( delimiters, tokenEnd );
  }

  if( tokenBegin < strEnd )
  {
    tokens.emplace_back( str.substr( tokenBegin, strEnd-tokenBegin ));
  }
  else if( !preTrimStr && str.find_first_of( delimiters, strEnd - 1 ) != string::npos )
  {
    tokens.emplace_back( "" );
  }

  return tokens;
}

/**
 * @brief Subdivide the string in substrings by whitespaces separators (see std::isspace()).
 *        Do not create any empty substrings.
 * @tparam CONTAINER The templated class of the results container (stdVector by default).
 * @param str The string to subdivide.
 * @return CONTAINER< string > The list of the subdivided substrings (stdVector< string > for instance).
 */
template< template< class ... > class CONTAINER = stdVector >
CONTAINER< string > tokenizeBySpaces( string const & str )
{
  return tokenize< CONTAINER >( str, " \f\n\r\t\v", true, true );
}

/**
 * @brief Trim the string
 * @param[in] str the string to trim
 * @param[in] charsToRemove the list of characters to remove
 * @return the trimmed string
 */
string_view trim( string_view str,
                  string_view charsToRemove );

/**
 * @brief Trim the left string
 * @param[in] s the string to trim
 * @return the trimmed string
 */
string_view ltrimSpaces( string_view s );

/**
 * @brief Trim the string so it does not starts nor ends with any whitespaces
 * @param[in] str the string to trim
 * @return the trimmed string
 */
string_view trimSpaces( string_view str );

/**
 * @brief Search for a string in the line, and return the line truncated before the string
 * @param[in] str the line to truncate
 * @param[in] strToRemove the string to search for in the line
 * @return the new (truncated) string
 */
string removeStringAndFollowingContent( string_view str,
                                        string_view strToRemove );

/**
 * @brief Add comma separators to an integral number for readability.
 * @tparam T the integral type of the number to format.
 * @param[in] num the integral number to format.
 * @return a string representation of the number with comma separators.
 */
template< typename T >
string addCommaSeparators( T const & num );

/**
 * @brief Divides a string by newline characters and returns a vector of strings containing each line.
 * Also calculates the width of the widest line.
 * @param linesWidth [out] Reference to a size_t that will be set to the width of the widest line
 * @param value The input string to divide into lines
 * @tparam STRING_T The type of the string (string or string_view)
 * @return A vector of STRING_T objects, each containing a single line from the input
 */
template< typename STRING_T >
stdVector< STRING_T > divideLines( size_t & linesWidth, string_view value );

/**
 * @brief Format all the lines by detecting spaces and by dividing each lines with maximum length specified.
 * If a word has a greater size than maxLength, it will be cut in one or many parts.
 * @param lines Vector containing all the lines to be formatted.
 * @param maxLineLength [inout] The max length a line can have.
 *                      The value is then set to the effective maximum line length
 * @tparam STRING_T The type of the string (string or string_view)
 * @return A vector containing the lines wrapped.
 */
template< typename STRING_T >
stdVector< STRING_T > wrapTextToMaxLength( stdVector< STRING_T > const & lines,
                                           size_t & maxLineLength );

/**
 * @brief Take a string, and return a array1d with the cast values
 * @tparam T the type to which the string will be cast
 * @param[in] str the string to turn into an array1d
 * @return the array1d that stores the cast values
 */
template< typename T >
array1d< T > fromStringToArray( string const & str )
{
  array1d< T > v;
  T sub;

  std::istringstream iss( str );
  while( iss >> sub )
  {
    v.emplace_back( sub );
  }
  return v;
}

/**
 * @brief Take a numerical value and convert/scale it to a string with a metric
 *  prefix. i.e. Kilo, Mega, Giga, Tera, Peta, Exa
 *
 * @tparam T Type of the value to be converted
 * @param value The value to be converted
 * @return String containging the scaled value.
 */
template< typename T >
string toMetricPrefixString( T const & value );

/**
 * @return The length of a constant string computed at compile-time.
 * @param str The null-character terminated constant string
 * @todo c++17: this function is to remove in favor of std::string_view
 */
constexpr size_t cstrlen( char const * const str )
{
  if( str )
  {
    char const * ptr = str;
    for(; *ptr != '\0'; ++ptr )
    {}
    return ptr - str;
  }
  else
  {
    return 0;
  }
}

/**
 * @return true if the string starts with the prefix.
 * @param str The string to compare
 * @param prefix A prefix we want to know if the string starts with.
 */
constexpr bool startsWith( std::string_view str, std::string_view prefix )
{
  return str.size() >= prefix.size() &&
         str.compare( 0, prefix.size(), prefix ) == 0;
}

/**
 * @return true if the string ends with the suffix.
 * @param str The string to compare
 * @param suffix A suffix we want to know if the string ends with.
 */
constexpr bool endsWith( std::string_view str, std::string_view suffix )
{
  return str.size() >= suffix.size() &&
         str.compare( str.size()-suffix.size(), suffix.size(), suffix ) == 0;
}

/**
 * @brief Overloading operator (<<) for std::optional<T>.
 *
 * This function displays the value contained in a std::optional<T> object if one exists.
 * Otherwise, it produces no output.
 *
 * @tparam T The type of the value contained std::optional.
 * @param os An output stream (for example, std::cout).
 * @param optValue std::optional<T> value to display.
 * @return The output stream
 */
template< typename T >
std::ostream & operator<<( std::ostream & os, std::optional< T > const & optValue )
{
  if( optValue )
  {
    os << optValue.value();
  }
  return os;
}

} // namespace stringutilities
} // namespace geos

#endif /* GEOS_COMMON_FORMAT_STRINGUTILITIES_HPP_ */
