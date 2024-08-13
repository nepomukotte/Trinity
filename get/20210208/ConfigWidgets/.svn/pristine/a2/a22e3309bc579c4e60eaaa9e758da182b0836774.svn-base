#include "ConfigTreeColumn.h"
#include "ConfigTreeRow.h"
#include "../UndoableConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "CCfg/Attribute.h"

namespace CCfg
{
namespace Model
{
namespace Column
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Column::ConfigTreeColumn class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param accessor the object for accessing the configuration model.
 */
ConfigTreeColumn::ConfigTreeColumn(UndoableConfigModelAccessor* accessor)
	: accessor(accessor)
{
}

/**
 * Destructor.
 */
ConfigTreeColumn::~ConfigTreeColumn()
{
}

/*
 *
 */
QColor ConfigTreeColumn::getRowForegroundTextColor(const AbstractRow& row) const
{
	const ConfigTreeRow& configRow = static_cast<const ConfigTreeRow&>(row);
	if(configRow.node->getBacktrackLevel() == 0)
	{
		if(configRow.node->isOverride())
		{
			return Qt::darkRed;
		}
		else
		{
			return Qt::black;
		}
	}
	else
	{
		return Qt::blue;
	}
}

/*
 *
 */
QColor ConfigTreeColumn::getRowBackgroundColor(const AbstractRow& row) const
{
	const ConfigTreeRow& configRow = static_cast<const ConfigTreeRow&>(row);
	if(configRow.node->hasValidValue())
	{
		return AbstractColumn::getRowBackgroundColor(row);
	}
	else
	{
		return QColor(0xFF, 0x80, 0x80);
	}
}

/*
 *
 */
Qt::ItemFlags ConfigTreeColumn::getRowFlags(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	if(node->isObject())
	{
		return Qt::ItemIsDropEnabled | Qt::ItemIsEnabled;
	}
	else
	{
		return Qt::ItemIsEnabled;
	}
}

}
}
}
