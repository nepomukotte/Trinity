#include "AbstractRow.h"

namespace CCfg
{
namespace Model
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::AbstractRow class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Destructor.
 */
AbstractRow::~AbstractRow()
{
}

/**
 * Indicates whether this row is editable.
 * It does not presume of the editability of a particular column, but only
 * indicates if something at the row level authorise or oppose to edition.
 * @return Returns true if the row may be edited.
 */
bool AbstractRow::isEditable() const
{
	return false;
}

}
}
