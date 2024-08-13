#ifndef CCFG_MODEL_COLUMN_CONFIGTREEROW_H_INCLUDED
#define CCFG_MODEL_COLUMN_CONFIGTREEROW_H_INCLUDED

#include "../AbstractRow.h"

namespace CCfg
{
namespace Model
{
// forward declaration
class ConfigNode;

namespace Column
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * Carries data about a row of the configuration tree model.
 *
 * @author fchateau
 */
struct ConfigTreeRow : public AbstractRow
{
	ConfigTreeRow(ConfigNode* node);
	virtual ~ConfigTreeRow();

	ConfigNode* node;

	virtual bool isEditable() const;
};

}
}
}

#endif // CCFG_MODEL_COLUMN_CONFIGTREEROW_H_INCLUDED
