#include "AttributeInclude.h"
#include "ConfigTreeRow.h"
#include "../UndoableConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "CCfg/Io/Url.h"

namespace CCfg
{
namespace Model
{
namespace Column
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Column::AttributeInclude class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
AttributeInclude::AttributeInclude(UndoableConfigModelAccessor* accessor)
	: ConfigTreeColumn(accessor)
{
}

/**
 * Destructor.
 */
AttributeInclude::~AttributeInclude()
{
}

/*
 *
 */
QString AttributeInclude::getHeaderDisplayName() const
{
	return "Include path";
}

/*
 *
 */
QString AttributeInclude::getRowDisplayData(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	if(node->isObject())
	{
		// Note: supporting general URL would require to do this only if the URL is a file
		return Io::Url::PathNormalizer::unescape(node->getIncludeUrl()).c_str();
	}
	else
	{
		return QString();
	}
}

/*
 *
 */
Qt::ItemFlags AttributeInclude::getRowFlags(const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	Qt::ItemFlags flags = ConfigTreeColumn::getRowFlags(row);
	if(row.isEditable() and node->isObject() and not node->isRoot())
	{
		flags |= Qt::ItemIsEditable;
	}
	return flags;
}

/*
 *
 */
bool AttributeInclude::setRowEditionData(const QString& newData, const AbstractRow& row) const
{
	ConfigNode* node = static_cast<const ConfigTreeRow&>(row).node;
	return accessor->setIncludeUrl(node, Io::Url::PathNormalizer::escape(newData.toStdString()));
}

}
}
}
