#include "AttributeValue.h"
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
// CCfg::Model::Column::AttributeValue class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param accessor the object for accessing the configuration model.
 */
AttributeValue::AttributeValue(UndoableConfigModelAccessor* accessor)
	: ConfigTreeColumn(accessor)
{
}

/**
 * Destructor.
 */
AttributeValue::~AttributeValue()
{
}

/*
 *
 */
QString AttributeValue::getHeaderDisplayName() const
{
	return "Value";
}

/*
 *
 */
QString AttributeValue::getRowDisplayData(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	if(not node->isObject())
	{
		return QString::fromStdString(node->getValue().getAsString());
	}
	else
	{
		return QString();
	}
}

/*
 *
 */
Qt::Alignment AttributeValue::getRowAlignment(const AbstractRow& /*row*/) const
{
	return Qt::AlignRight;
}

/*
 *
 */
Qt::ItemFlags AttributeValue::getRowFlags(const AbstractRow& row) const
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
bool AttributeValue::setRowEditionData(const QString& newData, const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	Value value;
	switch(node->getValue().getType())
	{
	case INTEGER:value.setInt(newData.toStdString());
		break;
	case HEXADECIMAL:value.setHex(newData.toStdString());
		break;
	case REAL:value.setReal(newData.toStdString());
		break;
	case BOOLEAN:value.setBool(newData.toStdString());
		break;
	case STRING:value.setString(newData.toStdString());
		break;
	default:
		break;
	}
	return accessor->setValue(node, value);
}

}
}
}
