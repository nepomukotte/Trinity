#ifndef CCFG_MODEL_COLUMN_ATTRIBUTEVALUE_H_INCLUDED
#define CCFG_MODEL_COLUMN_ATTRIBUTEVALUE_H_INCLUDED

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
 * Manages access to the "Value" column in the configuration tree model.
 *
 * @author fchateau
 */
class AttributeValue : public ConfigTreeColumn
{
public:
	AttributeValue(UndoableConfigModelAccessor* accessor);
	virtual ~AttributeValue();

	QString getHeaderDisplayName() const;
	QString getRowDisplayData(const AbstractRow& row) const;
	Qt::Alignment getRowAlignment(const AbstractRow& row) const;
	Qt::ItemFlags getRowFlags(const AbstractRow& row) const;

	bool setRowEditionData(const QString& newData, const AbstractRow& row) const;
};

}
}
}

#endif // CCFG_MODEL_COLUMN_ATTRIBUTEVALUE_H_INCLUDED
