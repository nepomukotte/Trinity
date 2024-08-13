#include "AbstractEntityBatch.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::AbstractEntityBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
AbstractEntityBatch::AbstractEntityBatch()
{
}

/**
 * Destructor.
 */
AbstractEntityBatch::~AbstractEntityBatch()
{
}

/**
 * @fn void AbstractEntityBatch::append(const AbstractEntity& entity)
 * Appends an entity on this batch.
 * The entity is copied, so it can be destroyed or overwritten immediately after
 * the call.
 */

/**
 * @fn size_t AbstractEntityBatch::size() const
 * Returns the number of objects in the batch.
 */

/**
 * @fn void AbstractEntityBatch::resize(size_t size)
 * Resizes the batch arrays to be ready for data incomming transfer.
 */

/**
 * @fn UidBatch& AbstractEntityBatch::getUidBatch(int columnIdx)
 * Returns the data of the specified column as a vector of UIDs.
 * @throw std::logic_error if the column is not of type UID. This should never
 * happen because the column index has no reason to be anything but a constant.
 * This exception is not meant to be catched, but to close the program like an
 * assertion.
 */

/**
 * @fn IntBatch& AbstractEntityBatch::getIntBatch(int columnIdx)
 * Returns the data of the specified column as a vector of integers.
 * @throw std::logic_error if the column is not of type string.
 */

/**
 * @fn StringBatch& AbstractEntityBatch::getStringBatch(int columnIdx)
 * Returns the data of the specified column as a vector of strings.
 * @throw std::logic_error if the column is not of type string.
 */

/**
 * @fn DateBatch& AbstractEntityBatch::getDateBatch(int columnIdx)
 * Returns the data of the specified column as a vector of dates.
 * @throw std::logic_error if the column is not of type date.
 */

}
}
}
