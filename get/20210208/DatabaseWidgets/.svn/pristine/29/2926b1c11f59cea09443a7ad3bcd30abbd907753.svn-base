/*
 * @file DatabaseFilterProxyModel.cpp
 * $Id$
 */

#include "DatabaseFilterProxyModel.h"
#include "AbstractDatabaseItem.h"
#include <cassert>

namespace CCfg
{
namespace Model
{
/**
 * Constructor.
 */
DatabaseFilterProxyModel::DatabaseFilterProxyModel(QObject* parent, FilterMode mode) : QSortFilterProxyModel(parent), filterMode(NoFilter)
{
	setFilterMode(mode);
}

/**
 * Sets filter mode.
 */
void DatabaseFilterProxyModel::setFilterMode(FilterMode newMode)
{
	filterMode = newMode;
}

/**
 * Filters out either:
 *  - nothing;
 *  - configurations and their versions;
 *  - runs and their conditions.
 */
bool DatabaseFilterProxyModel::filterAcceptsRow(int sourceRow, const QModelIndex & sourceParentIndex) const
{
	if (not sourceParentIndex.isValid()) return true;
	QModelIndex sourceIndex = sourceModel()->index(sourceRow, 0, sourceParentIndex);
	AbstractDatabaseItem* item = reinterpret_cast<AbstractDatabaseItem*>(sourceIndex.internalPointer());
	assert(item != 0);

	switch(item->getType())
	{
	case AbstractDatabaseItem::Config:
		if (filterMode & RunOnly) return false;
		break;
	case AbstractDatabaseItem::Version:
		if (filterMode & RunOnly) return false;
		break;
	case AbstractDatabaseItem::Site:
		if (filterMode & ConfigOnly) return false;
		break;
	case AbstractDatabaseItem::Run:
		if (filterMode & ConfigOnly) return false;
		break;
	case AbstractDatabaseItem::RunCondition:
		if (filterMode & ConfigOnly) return false;
		break;
	default:
		return true;
	}
	return true;
}

}
}
