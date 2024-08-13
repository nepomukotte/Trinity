#ifndef CCFG_MODEL_COLUMN_ATTRIBUTEINCLUDE_H_INCLUDED
#define CCFG_MODEL_COLUMN_ATTRIBUTEINCLUDE_H_INCLUDED

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
 * Manages access to the "Include Url" column in the configuration tree model.
 *
 * @author fchateau
 */
class AttributeInclude : public ConfigTreeColumn
{
public:
	AttributeInclude(UndoableConfigModelAccessor* accessor);
	virtual ~AttributeInclude();

	QString getHeaderDisplayName() const;
	QString getRowDisplayData(const AbstractRow& row) const;
	Qt::ItemFlags getRowFlags(const AbstractRow& row) const;

	bool setRowEditionData(const QString& newData, const AbstractRow& row) const;
};
}
}
}

#endif // CCFG_MODEL_COLUMN_ATTRIBUTEINCLUDE_H_INCLUDED
