#include "AttributeName.h"
#include "ConfigTreeRow.h"
#include "../UndoableConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "CCfg/PathNode.h"
#include "CCfg/Comment.h"

namespace CCfg
{
namespace Model
{
namespace Column
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Column::AttributeName class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param accessor The object for accessing the configuration model.
 */
AttributeName::AttributeName(UndoableConfigModelAccessor* accessor)
	: ConfigTreeColumn(accessor)
{
}

/**
 * Destructor.
 */
AttributeName::~AttributeName()
{
}

/*
 *
 */
QString AttributeName::getHeaderDisplayName() const
{
	return "Name";
}

/*
 *
 */
QString AttributeName::getRowDisplayData(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	return QString::fromStdString(node->getFullName());
}

/*
 *
 */
QString AttributeName::getRowTooltip(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	std::string tooltip;
	tooltip += node->getPreComment().getText();
	if(not node->getPreComment().getText().empty() and
		not node->getPostComment().getText().empty())
	{
		tooltip += '\n';
	}
	tooltip += node->getPostComment().getText();
	return QString::fromStdString(tooltip);
}

/*
 *
 */
QFont AttributeName::getRowFont(const AbstractRow& row) const
{
	QFont font;
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	if(node->hasScalarInstances())
	{
		font.setBold(true);
		// no need to initialize the font completely, only boldness will be
		// applied on the default font
	}
	return font;
}

/*
 *
 */
Qt::ItemFlags AttributeName::getRowFlags(const AbstractRow& row) const
{
	Qt::ItemFlags flags = Qt::ItemIsSelectable
		| Qt::ItemIsDragEnabled
		| ConfigTreeColumn::getRowFlags(row);

	if(row.isEditable())
	{
		flags |= Qt::ItemIsEditable;
	}
	return flags;
}

/*
 *
 */
bool AttributeName::setRowEditionData(const QString& newData, const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	return accessor->setNameAndIndex(node, newData.toStdString());
}

}
}
}
