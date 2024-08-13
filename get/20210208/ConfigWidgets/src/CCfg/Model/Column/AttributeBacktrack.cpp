#include "AttributeBacktrack.h"
#include "ConfigTreeRow.h"
#include "../ConfigNode.h"

namespace CCfg
{
namespace Model
{
namespace Column
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Column::AttributeBacktrack class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
AttributeBacktrack::AttributeBacktrack()
	: ConfigTreeColumn(0)
{
}

/**
 * Destructor.
 */
AttributeBacktrack::~AttributeBacktrack()
{
}

/*
 *
 */
QString AttributeBacktrack::getHeaderDisplayName() const
{
	return "BL";
}

QString AttributeBacktrack::getHeaderTooltip() const
{
	return QObject::tr("<b>The Backtrack Level</b>"
	"<p>For an inherited attribute, it is the number of non-inherited parent "
	"objects that must be crossed before reaching the first common ancestor "
	"with the concrete attribute.</p>");
}

/*
 *
 */
QString AttributeBacktrack::getRowDisplayData(const AbstractRow& row) const
{
	int btLevel = static_cast<const ConfigTreeRow&>(row).node->getBacktrackLevel();
	if(btLevel != 0)
	{
		return QString::number(btLevel);
	}
	else
	{
		return QString();
	}
}

/*
 *
 */
Qt::Alignment AttributeBacktrack::getRowAlignment(const AbstractRow& /*row*/) const
{
	return Qt::AlignRight;
}

}
}
}
