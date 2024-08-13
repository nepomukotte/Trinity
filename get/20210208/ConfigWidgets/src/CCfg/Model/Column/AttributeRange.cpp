#include "AttributeRange.h"
#include "ConfigTreeRow.h"
#include "../UndoableConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "CCfg/Range.h"

namespace CCfg
{
namespace Model
{
namespace Column
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Column::AttributeRange class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param accessor the object for accessing the configuration model.
 */
AttributeRange::AttributeRange(UndoableConfigModelAccessor* accessor)
	: ConfigTreeColumn(accessor)
{
}

/**
 * Destructor.
 */
AttributeRange::~AttributeRange()
{
}

/*
 *
 */
QString AttributeRange::getHeaderDisplayName() const
{
	return "Range";
}

/*
 *
 */
QString AttributeRange::getRowDisplayData(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	if(not node->isObject())
	{
		Range* range = node->getRange();
		if(range != 0)
		{
			return QString::fromStdString(range->toString());
		}
	}
	return QString();
}

/*
 *
 */
Qt::ItemFlags AttributeRange::getRowFlags(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;

	Qt::ItemFlags flags = ConfigTreeColumn::getRowFlags(row);
	if(row.isEditable() and not node->isObject())
	{
		flags |= Qt::ItemIsEditable;
	}
	return flags;
}

/*
 *
 */
bool AttributeRange::setRowEditionData(const QString& newData, const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	return accessor->setRange(node, Range::createFromString(newData.toStdString()));
}

}
}
}
