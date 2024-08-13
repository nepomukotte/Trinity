#include "AbstractDatabaseItem.h"

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::AbstractDatabaseItem class methods body
////////////////////////////////////////////////////////////////////////////////

const size_t AbstractDatabaseItem::INVALID_INDEX;

/**
 * Default constructor.
 */
AbstractDatabaseItem::AbstractDatabaseItem()
{
}

/**
 * Destructor.
 */
AbstractDatabaseItem::~AbstractDatabaseItem()
{
}

/**
 * Searches the specified item among the children of this one.
 * @return Returns the index of the item, or
 */
size_t AbstractDatabaseItem::searchItem(AbstractDatabaseItem* item)
{
	for(size_t i = 0, n = getItemCount(); i < n; ++i)
	{
		if(getItem(i) == item)
		{
			return i;
		}
	}
	return INVALID_INDEX;
}

/**
 * Returns whether this item is a terminal item (cannot have children) or not.
 */
bool AbstractDatabaseItem::isLeaf() const
{
	return false; // use false, because most items will use this value.
}

/**
 * @fn bool AbstractDatabaseItem::canSetName() const
 * Returns whether it is possible to change the name of this item using setName.
 */

/**
 * @fn void AbstractDatabaseItem::setName(const std::string& newName)
 * Changes the name of this database item, if possible.
 * @throw std::runtime_error if the name cannot be be (ie: canSetName returns false).
 */

/**
 * @fn const std::string& AbstractDatabaseItem::getName() const
 * Returns the edit name of the database item.
 */

/**
 * Returns the name to display for this item.
 * The default implementation returns the general name.
 */
const std::string& AbstractDatabaseItem::getDisplayName() const
{
	return getName();
}

/**
 * @fn size_t AbstractDatabaseItem::getItemCount()
 * Returns the current number of children of this item.
 */

/**
 * @fn AbstractDatabaseItem* AbstractDatabaseItem::getItem(size_t index)
 * Returns the child item at the specified index.
 * @precondition The index must be between 0 and getItemCount()-1 included.
 */

/**
 * @fn AbstractDatabaseItem* AbstractDatabaseItem::getParent()
 * Returns the parent of this item.
 */

/**
 * @fn void AbstractDatabaseItem::refresh()
 * Refresh the list of children items.
 */

/**
 * @fn void AbstractDatabaseItem::removeItem(size_t index, RemoveOptions options)
 * Removes the child item at the specified index.
 * @param index the index of the children to remove.
 * @param options determines the extent of the removal.
 * If the flag FromList is set the item is removed from the children list.
 * If the flag FromDataStore is set, the removal is permanent: the underlying
 * representation is erased.
 * Theses flags are independant, the 3 combinations are allowed.
 */

/**
 * @fn void AbstractDatabaseItem::clear()
 * Clears all children items.
 */

/**
 * @fn std::string AbstractDatabaseItem::buildUrl() const
 * Builds a URL for this item.
 */

/**
 * @fn AbstractDatabaseItem::ItemType AbstractDatabaseItem::getType() const
 * Returns the type of the concrete item.
 */

}
}
