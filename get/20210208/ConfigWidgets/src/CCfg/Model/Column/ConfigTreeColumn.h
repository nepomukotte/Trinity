#ifndef CCFG_MODEL_COLUMN_CONFIGTREECOLUMN_H_INCLUDED
#define CCFG_MODEL_COLUMN_CONFIGTREECOLUMN_H_INCLUDED

#include "../AbstractColumn.h"

namespace CCfg
{
namespace Model
{
// forward declarations
class UndoableConfigModelAccessor;

namespace Column
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * Provides the shared methods implementation for classes that manages
 * columns of the configuration model.
 *
 * @author fchateau
 */
class ConfigTreeColumn : public AbstractColumn
{
protected:
	/// the object for accessing the configuration model.
	UndoableConfigModelAccessor* accessor;

public:
	ConfigTreeColumn(UndoableConfigModelAccessor* accessor);
	virtual ~ConfigTreeColumn();

	virtual QString getHeaderDisplayName() const = 0;
	virtual QColor getRowForegroundTextColor(const AbstractRow& row) const;
	virtual QColor getRowBackgroundColor(const AbstractRow& row) const;
	virtual Qt::ItemFlags getRowFlags(const AbstractRow& row) const;
};

}
}
}

#endif // CCFG_MODEL_COLUMN_CONFIGTREECOLUMN_H_INCLUDED
