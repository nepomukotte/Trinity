#ifndef CCFG_MODEL_COLUMN_ATTRIBUTEBACKTRACK_H_INCLUDED
#define CCFG_MODEL_COLUMN_ATTRIBUTEBACKTRACK_H_INCLUDED

#include "ConfigTreeColumn.h"

namespace CCfg
{
namespace Model
{
namespace Column
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * Manages access to the "backtrack level" column in the configuration tree model.
 *
 * @author fchateau
 */
class AttributeBacktrack : public ConfigTreeColumn
{
public:
	AttributeBacktrack();
	virtual ~AttributeBacktrack();

	QString getHeaderDisplayName() const;
	QString getHeaderTooltip() const;
	QString getRowDisplayData(const AbstractRow& row) const;
	Qt::Alignment getRowAlignment(const AbstractRow& row) const;
};

}
}
}

#endif // CCFG_MODEL_COLUMN_ATTRIBUTEBACKTRACK_H_INCLUDED
