#include "AttributeUnit.h"
#include "ConfigTreeRow.h"
#include "../UndoableConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "CCfg/Value.h"

namespace CCfg
{
namespace Model
{
namespace Column
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Column::AttributeUnit class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param accessor the object for accessing the configuration model.
 */
AttributeUnit::AttributeUnit(UndoableConfigModelAccessor* accessor)
	: ConfigTreeColumn(accessor)
{
}

/**
 * Destructor.
 */
AttributeUnit::~AttributeUnit()
{
}

/*
 *
 */
QString AttributeUnit::getHeaderDisplayName() const
{
	return "Unit";
}

/*
 *
 */
QString AttributeUnit::getRowDisplayData(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	if(not node->isObject())
	{
		return QString::fromStdString(node->getUnit());
	}
	else
	{
		return QString();
	}
}

/*
 *
 */
Qt::ItemFlags AttributeUnit::getRowFlags(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;

	Qt::ItemFlags flags = ConfigTreeColumn::getRowFlags(row);
	if(row.isEditable() and not node->isObject())
	{
		ValueTypeId type = node->getValue().getType();
		if(type == INTEGER or type == HEXADECIMAL or type == REAL)
		{
			flags |= Qt::ItemIsEditable;
		}
	}
	return flags;
}

/*
 *
 */
bool AttributeUnit::setRowEditionData(const QString& newData, const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	return accessor->setUnit(node, newData.toStdString());
}

}
}
}
