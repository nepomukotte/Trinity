#include "AbstractColumn.h"
#include "AbstractRow.h"

namespace CCfg
{
namespace Model
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::AbstractColumn class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Destructor.
 */
AbstractColumn::~AbstractColumn()
{
}

/**
 * @fn QString AbstractColumn::getHeaderDisplayName() const
 * Provides the display name of the column header.
 * @return Returns the display name of the header of this column.
 */

/**
 * Provides the tooltip text for the column header.
 * @return Returns the tooltip text for the column header.
 */
QString AbstractColumn::getHeaderTooltip() const
{
	return QString();
}

/**
 * @fn QString AbstractColumn::getRowDisplayData(const AbstractRow& row) const
 * Provides the display text of the data for the specified row.
 * @param row the row this method applies on.
 * @return Returns the display text of the data for the specified row.
 */

/**
 * Provides the edition data for the specified row.
 * @param row the row this method applies on.
 * @return Returns the edition data for the specified row.
 */
QString AbstractColumn::getRowEditionData(const AbstractRow& row) const
{
	return getRowDisplayData(row);
}

/**
 * Provides the tooltip text for the specified row.
 * @param row the row this method applies on.
 * @return Returns the tooltip text for the specified row.
 */
QString AbstractColumn::getRowTooltip(const AbstractRow& /*row*/) const
{
	return QString();
}

/**
 * Provides the font for the text of the specified row.
 */
QFont AbstractColumn::getRowFont(const AbstractRow& /*row*/) const
{
	return QFont();
}

/**
 * Provides the data alignment for the specified row.
 * @param row the row this method applies on.
 * @return Returns the data alignment for the specified row.
 */
Qt::Alignment AbstractColumn::getRowAlignment(const AbstractRow& /*row*/) const
{
	return Qt::AlignLeft;
}

/**
 * Provides the foreground text color for the specified row.
 * @param row the row this method applies on.
 * @return Returns the foreground text color for the specified row.
 */
QColor AbstractColumn::getRowForegroundTextColor(const AbstractRow& /*row*/) const
{
	return Qt::black;
}

/**
 * Provides the background color for the specified row.
 * @param row the row this method applies on.
 * @return Returns the background color for the specified row.
 */
QColor AbstractColumn::getRowBackgroundColor(const AbstractRow& /*row*/) const
{
	return Qt::white;
}

/**
 * Provides the flags for the specified row.
 * @param row the row this method applies on.
 */
Qt::ItemFlags AbstractColumn::getRowFlags(const AbstractRow& /*row*/) const
{
	return Qt::ItemFlags();
}

/**
 * Defines the edition data for the specified row.
 * @param newData the new edition data.
 * @param row the row this method applies on.
 * @return Returns true if edition was accepted, false otherwise.
 */
bool AbstractColumn::setRowEditionData(const QString& /*newData*/,
	const AbstractRow& /*row*/) const
{
	return false;
}


}
}
