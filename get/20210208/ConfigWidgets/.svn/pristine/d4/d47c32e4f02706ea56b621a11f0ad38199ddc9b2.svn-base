#include "ConfigTreeRow.h"
#include "../ConfigNode.h"
#include <assert.h>

namespace CCfg
{
namespace Model
{
namespace Column
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Column::ConfigTreeRow class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConfigTreeRow::ConfigTreeRow(ConfigNode* node)
	: node(node)
{
	assert(node != 0);
}

/**
 * Destructor.
 */
ConfigTreeRow::~ConfigTreeRow()
{
}

/*
 *
 */
bool ConfigTreeRow::isEditable() const
{
	return not node->isInherited();
}

}
}
}
