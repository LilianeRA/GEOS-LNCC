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
 * @file TableData.hpp
 */

#ifndef GEOS_COMMON_FORMAT_TABLE_TABLEDATA_HPP
#define GEOS_COMMON_FORMAT_TABLE_TABLEDATA_HPP

#include "common/Units.hpp"
#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"
#include "TableTypes.hpp"

namespace geos
{

/**
 * @brief Class for managing table data
 */
class TableData
{
public:

  /**
   * @brief Representing a data in TableData
   */
  struct CellData
  {
    /// The cell type
    CellType type;
    /// The cell value
    string value;

    /**
     * @brief Comparison operator for cell value
     * @param other The cell data value to compare
     * @return The comparison result
     */
    bool operator==( CellData const & other ) const
    {
      return value == other.value;
    }
  };

  /// Alias for table data rows with cells values
  using DataRows = stdVector< stdVector< CellData > >;

  /**
   * @brief Add a row to the table.
   * The values passed to addRow (can be any type).
   * @param args CellData values to be added to the row.
   */
  template< typename ... Args >
  void addRow( Args const & ... args );

  /**
   * @brief Add a row to the table
   * @param row A vector of string representing a row
   */
  void addRow( stdVector< CellData > const & row );

  /**
   * @brief Add a line separator to the table
   * You must have filled values in TableData before using it
   */
  void addSeparator();

  /**
   * @brief Reset data in the table
   */
  void clear();

  /**
   * @return The rows of the table
   */
  stdVector< stdVector< CellData > > const & getTableDataRows() const;

  /**
   * @brief Get all error messages
   * @return The vector of error messages
   */
  stdVector< string > const & getErrorMsgs() const;

  /**
   * @return The const table data rows
   */
  DataRows const & getCellsData() const
  { return m_rows; }

  /**
   * @brief Comparison operator for data rows
   * @param comparingTable The tableData values to compare
   * @return The comparison result
   */
  inline bool operator==( TableData const & comparingTable ) const
  {

    return getCellsData() == comparingTable.getCellsData();
  }

private:
  /// vector containing all rows with cell values
  DataRows m_rows;

};


/**
 * @brief Class for managing 2D table m_data
 */
class TableData2D
{
public:

  /// Type real64 for a row
  using RowType = real64;
  /// Type real64 for a column
  using ColumnType = real64;

  /// Struct containing conversion informations
  struct TableDataHolder
  {
    /// Vector containing all columns names
    /// A header value is presented as "pressure [K] = {}"
    stdVector< string > headerNames;
    /// TableData to be built
    TableData tableData;
  };

  /**
   * @brief Add a cell to the table. If necessary, create automatically the containing column & row.
   * @tparam T The value passed to addCell (can be any type).
   * @param value CellData value to be added.
   * @param rowValue The value of the row containing the cell.
   * @param columnValue The value of the column containing the cell.
   */
  template< typename T >
  void addCell( RowType rowValue, ColumnType columnValue, T const & value );

  /**
   * @brief Collects all the values needed to build the table
   * @param dim0AxisCoordinates Vector containing all row axis values
   * @param dim1AxisCoordinates Vector containing all column axis values
   * @param values Array containing all table values contiguously
   * @param columnMajorValues Set the row/column major convention
   */
  void collectTableValues( arraySlice1d< real64 const > dim0AxisCoordinates,
                           arraySlice1d< real64 const > dim1AxisCoordinates,
                           arrayView1d< real64 const > values,
                           bool columnMajorValues );

  /**
   * @param values Vector containing all table values
   * @param valueDescription The description of the value (typically, the value unit description)
   * @param columnMajorValues Set the row/column major convention
   * @param coordinates Array containing row/column axis values
   * @param rowAxisDescription The description for a row unit value
   * @param columnAxisDescription The description for a column unit value
   * @return A struct containing the tableData converted and all header values ;
   */
  TableData2D::TableDataHolder convertTable2D( ArrayOfArraysView< real64 const > const coordinates,
                                               string_view rowAxisDescription,
                                               string_view columnAxisDescription,
                                               arrayView1d< real64 const > const values,
                                               bool columnMajorValues,
                                               string_view valueDescription );

  /**
   * @return Convert and return a struct containing a 1D Table, the column names list from a TableData2D and any errors related to the table
   * @param dataDescription The table dataDescription shown at the top left side
   * @param rowFmt The y axis units of the table.
   * @param columnFmt  The x axis units of the table.
   * @note The rows and columns FMT can be customized. The bracket "{}" will be replaced by the axis value.
   * By default it displays the axis value.
   * I.E to display a customized axis to show the pressures in y axis, a rowFmt value can be : "pressure [K] = {}"
   */
  TableDataHolder buildTableData( string_view dataDescription,
                                  string_view rowFmt = "{}", string_view columnFmt = "{}" ) const;

private:
  /// @brief all cell values by their [ row ][ column ]
  std::map< RowType, std::map< ColumnType, string > > m_data;

  /// @brief Store all column values when adding cell
  std::set< real64 > m_columnValues;
};

/**
 * @brief Trait to check is the args is a special type of cell
 * @tparam T The type of a cell
 */
template< typename T >
constexpr bool isCellType = std::is_same_v< T, CellType >;

template< typename ... Args >
void TableData::addRow( Args const &... args )
{
  stdVector< CellData > cells;
  ( [&] {
    static_assert( has_formatter_v< decltype(args) > || isCellType< std::decay_t< decltype(args) > >, "Argument passed in addRow cannot be converted to string nor a CellType" );
    if constexpr (std::is_same_v< Args, CellType >) {
      cells.push_back( { args, string() } );
    }
    else
    {
      cells.push_back( {CellType::Value, GEOS_FMT( "{}", args )} );
    }
  } (), ...);
  addRow( cells );
}

template< typename T >
void TableData2D::addCell( real64 const rowValue, real64 const columnValue, T const & value )
{
  static_assert( has_formatter_v< decltype(value) >, "Argument passed in addCell cannot be converted to string" );
  m_columnValues.insert( columnValue );
  m_data[rowValue][columnValue] = GEOS_FMT( "{}", value );
}

}
#endif /* GEOS_COMMON_FORMAT_TABLE_TABLEDATA_HPP */
