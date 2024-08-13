#ifndef CCFG_MODEL_ABSTRACTROW_H_INCLUDED
#define CCFG_MODEL_ABSTRACTROW_H_INCLUDED

namespace CCfg
{
namespace Model
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This is the base class describing rows of a model.
 * By deriving it, one can pass arbitrary data to the methods of AbstractColumn.
 *
 * @author fchateau
 */
struct AbstractRow
{
	virtual ~AbstractRow();

	virtual bool isEditable() const;
};

}
}

#endif // CCFG_MODEL_ABSTRACTROW_H_INCLUDED
