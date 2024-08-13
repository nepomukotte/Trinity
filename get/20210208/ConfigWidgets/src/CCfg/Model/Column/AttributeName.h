#ifndef CCFG_MODEL_COLUMN_ATTRIBUTENAME_H_INCLUDED
#define CCFG_MODEL_COLUMN_ATTRIBUTENAME_H_INCLUDED

#include "ConfigTreeColumn.h"

namespace CCfg
{
namespace Model
{
namespace Column
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * Manages access to the "Name" column in the configuration tree model.
 *
 * @author fchateau
 */
class AttributeName : public ConfigTreeColumn
{
public:
	AttributeName(UndoableConfigModelAccessor* accessor);
	virtual ~AttributeName();

	QString getHeaderDisplayName() const;
	QString getRowDisplayData(const AbstractRow& row) const;
	QString getRowTooltip(const AbstractRow& row) const;
	QFont getRowFont(const AbstractRow& row) const;
	Qt::ItemFlags getRowFlags(const AbstractRow& row) const;

	bool setRowEditionData(const QString& newData, const AbstractRow& row) const;
};

}
}
}

#endif // CCFG_MODEL_COLUMN_ATTRIBUTENAME_H_INCLUDED
