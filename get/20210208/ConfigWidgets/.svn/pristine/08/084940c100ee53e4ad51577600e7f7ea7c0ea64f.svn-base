#ifndef CCFG_MODEL_COLUMN_ATTRIBUTERANGE_H_INCLUDED
#define CCFG_MODEL_COLUMN_ATTRIBUTERANGE_H_INCLUDED

#include "ConfigTreeColumn.h"

namespace CCfg
{
namespace Model
{
namespace Column
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * Manages access to the "Range" column in the configuration tree model.
 *
 * @author fchateau
 */
class AttributeRange : public ConfigTreeColumn
{
public:
	AttributeRange(UndoableConfigModelAccessor* accessor);
	virtual ~AttributeRange();

	QString getHeaderDisplayName() const;
	QString getRowDisplayData(const AbstractRow& row) const;
	Qt::ItemFlags getRowFlags(const AbstractRow& row) const;

	bool setRowEditionData(const QString& newData, const AbstractRow& row) const;

private:
	AttributeRange(const AttributeRange& rvalue);
	AttributeRange& operator=(const AttributeRange& rvalue);
};

}
}
}

#endif // CCFG_MODEL_COLUMN_ATTRIBUTERANGE_H_INCLUDED
