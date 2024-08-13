#include "ConfigEditorTree.h"
#include "CommonParentSelectionModel.h"
#include "../Model/ConfigModel.h"
#include "../Model/ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include <QUndoStack>
#include <QMenu>
#include <vector>

namespace CCfg
{
using Model::ConfigModel;
using Model::ConfigNode;
using Model::IndexList;

namespace Ui
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::ConfigEditorTree class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param parent the parent widget.
 */
ConfigEditorTree::ConfigEditorTree(QWidget* parent)
	: QTreeView(parent)
{
	setAllColumnsShowFocus(true);
	setSelectionBehavior(QAbstractItemView::SelectRows);
	setSelectionMode(QAbstractItemView::ExtendedSelection);
}

/**
 * Destructor.
 */
ConfigEditorTree::~ConfigEditorTree()
{
}

/**
 * Sets the model for the view to present.
 */
void ConfigEditorTree::setModel(QAbstractItemModel* model)
{
	try
	{
		if(model == this->model())
		{
			return;
		}

		if(this->model() != 0)
		{
			disconnect(this->model(), 0, this, 0);
		}
		QTreeView::setModel(model);
		setSelectionModel(new CommonParentSelectionModel(model, this));

		setColumnWidth(ConfigModel::NameColumn, columnWidth(ConfigModel::NameColumn) * 2);
		resizeColumnToContents(ConfigModel::BacktrackColumn);

		connect(model, SIGNAL(beginUpdate()),
			this, SLOT(onBeginUpdate()));
		connect(model, SIGNAL(endUpdate()),
			this, SLOT(onEndUpdate()));
		connect(model, SIGNAL(endInsertion(ConfigNode*, IndexList)),
			this, SLOT(onEndInsertion(ConfigNode*, IndexList)));
		connect(model, SIGNAL(beginRemoval(ConfigNode*, IndexList)),
			this, SLOT(onBeginRemoval(ConfigNode*, IndexList)));
	}
	CATCH_AND_REPORT
}

/**
 * Returns the index of the configuration root.
 */
QModelIndex ConfigEditorTree::getConfigRoot() const
{
	return model()->index(0, 0, QModelIndex());
}

/**
 * Recursively expands all item starting from the specified item.
 */
void ConfigEditorTree::expandAll(const QModelIndex& item)
{
	QAbstractItemModel* mdl = model();
	expand(item);
	if(mdl->hasChildren(item))
	{
		int rowCount = mdl->rowCount(item);
		for(int i = 0; i < rowCount; ++i)
		{
			expandAll(mdl->index(i, 0, item));
		}
	}
}

/**
 * Recursively collapses all item starting from the specified item.
 */
void ConfigEditorTree::collapseAll(const QModelIndex& item)
{
	QAbstractItemModel* mdl = model();
	collapse(item);
	if(mdl->hasChildren(item))
	{
		int rowCount = mdl->rowCount(item);
		for(int i = 0; i < rowCount; ++i)
		{
			collapseAll(mdl->index(i, 0, item));
		}
	}
}

/**
 * Prepares to row removal by backing up their expand state.
 * @param parent the parent of the removed rows.
 * @param start the index of the first removed item.
 * @param end the index of the last removed item.
 */
void ConfigEditorTree::rowsAboutToBeRemoved(const QModelIndex& parent, int start, int end)
{
	try
	{
		saveSelection(parent);
		saveExpandingState(parent);
		QTreeView::rowsAboutToBeRemoved(parent, start, end);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to row insertion by restoring backed up expand state.
 * @param parent the parent of the inserted rows.
 * @param start the index of the first inserted item.
 * @param end the index of the last inserted item.
 */
void ConfigEditorTree::rowsInserted(const QModelIndex& parent, int start, int end)
{
	try
	{
		QTreeView::rowsInserted(parent, start, end);
		restoreExpandingState(parent);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the beginning of a composite update on the model.
 */
void ConfigEditorTree::onBeginUpdate()
{
	try
	{
		expandState.clear();
		selectionBackup.clear();
		currentSelectionBackup.clear();
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the end of a composite update on the model.
 */
void ConfigEditorTree::onEndUpdate()
{
	try
	{
		restoreSelection();
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to an insertion in the model.
 * @param parent the parent node of the newly inserted nodes.
 * @param indexes the indexes (relative to the parent) of the inserted nodes.
 */
void ConfigEditorTree::onEndInsertion(ConfigNode* parent, const IndexList& indexes)
{
	try
	{
		ConfigModel* mdl = static_cast<ConfigModel*>(model());
		QItemSelectionModel* selection = selectionModel();
		std::vector<ConfigNode>& children = parent->getChildren();

		selection->clearSelection();
		for(IndexList::const_iterator index = indexes.begin(); index != indexes.end(); ++index)
		{
			selection->select(mdl->nodeToIndex(&children[*index]),
				QItemSelectionModel::Select | QItemSelectionModel::Rows);
		}
		if(indexes.size() > 0)
		{
			selection->setCurrentIndex(mdl->nodeToIndex(&children[indexes.front()]),
				QItemSelectionModel::Current | QItemSelectionModel::Rows);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a removal in the model.
 * @param parent the parent node of the removed nodes.
 * @param indexes the indexes (relative to the parent) of the nodes that will be
 * removed.
 */
void ConfigEditorTree::onBeginRemoval(ConfigNode* parent, const IndexList& indexes)
{
	try
	{
		if(indexes.size() == 0)
		{
			return;
		}

		std::vector<ConfigNode>& children = parent->getChildren();
		ConfigNode* nodeToSelect;// the node to select

		// if all nodes have been removed
		if(children.size() == indexes.size())
		{
			// select the parent
			nodeToSelect = parent;
		}
		else
		{
			// otherwise, there are nodes left. They can be at the following 3 places
			// (in order of preference for selection):
			// - after the removed nodes;
			// - interleaved with the removed nodes;
			// - before the removed nodes;

			// if the last removed node is not the last children.
			// It means there is a least 1 nodes after.
			if(children.size() > indexes.back() + 1)
			{
				// select the node immediately after the last removed node
				nodeToSelect = &children[indexes.back() + 1];
			}
			// if we removed only 1 node, there can't be interleaved nodes, so it
			// must be placed before.
			else if(indexes.size() == 1)
			{
				assert(indexes.front() != 0);
				// select the node before the first removed node
				nodeToSelect = &children[indexes.front() - 1];
			}
			else
			{
				// look for interleaved nodes.

				IndexList::const_reverse_iterator curIndex = indexes.rbegin();
				IndexList::const_reverse_iterator prevIndex = ++indexes.rbegin();
				// search the last gap in the index list of the removed nodes.
				while(prevIndex != indexes.rend() and *curIndex == *prevIndex + 1)
				{
					++curIndex;
					++prevIndex;
				}
				// if the loop ended because we could not find a gap.
				if(prevIndex == indexes.rend())
				{
					assert(indexes.front() != 0);
					// select the node before the first removed node
					nodeToSelect = &children[indexes.front() - 1];
				}
				else // a gap was found.
				{
					// select the last node of the gap
					nodeToSelect = &children[*curIndex - 1];
				}
			}
		}
		ConfigModel* mdl = static_cast<ConfigModel*>(model());
		QModelIndex index = mdl->nodeToIndex(nodeToSelect);
		selectionModel()->setCurrentIndex(index, QItemSelectionModel::ClearAndSelect
			| QItemSelectionModel::Current | QItemSelectionModel::Rows);
	}
	CATCH_AND_REPORT
}


/**
 * Recursively backs up the path of expanded items under the specified node.
 * More precisely it backs up all children of the specified item.
 */
void ConfigEditorTree::saveExpandingState(const QModelIndex& parentIdx)
{
	ConfigModel* model = static_cast<ConfigModel*>(this->model());

	if(not parentIdx.isValid())
	{
		saveExpandingState(model->index(0, 0, parentIdx));
		return;
	}

	ConfigNode* node = model->indexToNode(parentIdx);

	if(isExpanded(parentIdx))
	{
		expandState.insert(node->getPath().toString());
	}

	if(node->isObject() and node->isPopulated())
	{
		std::vector<ConfigNode>& children = node->getChildren();
		for(size_t row = 0; row < children.size(); ++row)
		{
			QModelIndex child = model->index(row, ConfigModel::NameColumn, parentIdx);
			saveExpandingState(child);
		}
	}
}

/**
 * Expands the specified node and all its children being in the backup.
 * @param parentIdx the index of the parent of the items to restore.
 */
void ConfigEditorTree::restoreExpandingState(const QModelIndex& parentIdx)
{
	ConfigModel* model = static_cast<ConfigModel*>(this->model());
	std::set<std::string>::const_iterator begin;
	std::set<std::string>::const_iterator end;

	if(not parentIdx.isValid())
	{
		// defines the range in the set
		begin = expandState.begin();
		end = expandState.end();
	}
	else
	{
		ConfigNode* node = model->indexToNode(parentIdx);

		// We initialize 2 strings that define a range for the paths of all children nodes.
		// Example: for a parent path: A[0].B[1], the range is A[0].B[1^
		// ('^' is the character after ']').
		// Every paths beginning with A[0].B[1] fall within this range.
		// This is like searching a word in a dictionnary: if you want all words
		// beginning with "ar" you must use the range ["ar", "as"[ because longer
		// words like "armor", "art"... are lesser that "as" in lexicographical order.
		// This is why we increment the last character by 1.
		std::string lowerPath = node->getPath().toString();
		std::string upperPath = lowerPath;
		++upperPath[upperPath.length() - 1]; // increment the last character to close the range.

		// defines the range in the set
		begin = expandState.lower_bound(lowerPath);
		end = expandState.upper_bound(upperPath);
	}

	// browse the range to expand items.
	for(std::set<std::string>::const_iterator item = begin; item != end; ++item)
	{
		QModelIndex index = model->pathToIndex(Path(*item));
		setExpanded(index, true);
	}
}

/**
 * Saves the path of selected and current items in internal attributes.
 * @param parentIdx only selected items that are descendants of this parent are saved.
 */
void ConfigEditorTree::saveSelection(const QModelIndex& parentIdx)
{
	ConfigModel* mdl = static_cast<ConfigModel*>(model());

	// we cannot use selection ranges because selection boundaries can be deleted or renamed.
	QModelIndexList selection = selectionModel()->selectedRows();

	// only save selected items that are descendants of the specified parent.
	// Testing the first selected item is enough, because our selection model
	// does not allow selection from different parents.
	if(not selection.empty() and mdl->isDescendant(selection.front(), parentIdx))
	{
		foreach(QModelIndex index, selection)
		{
			selectionBackup.push_back(mdl->indexToNode(index)->getPath());
		}
		QModelIndex currentIndex = selectionModel()->currentIndex();
		if(currentIndex.isValid())
		{
			currentSelectionBackup = mdl->indexToNode(currentIndex)->getPath();
		}
	}
}

/**
 * Restores every selected item that were saved by previous calls to saveSelection.
 * If some selected items have disappeared, they will be ignored. If the saved
 * "current item" has disappeared, its parent (or the first existing ascendant)
 * will be used.
 */
void ConfigEditorTree::restoreSelection()
{
	ConfigModel* mdl = static_cast<ConfigModel*>(model());
	QItemSelectionModel* selection = selectionModel();

	selection->clearSelection();
	// selects every saved item path.
	for(std::list<Path>::iterator path = selectionBackup.begin(); path != selectionBackup.end(); ++path)
	{
		selection->select(mdl->pathToIndex(*path),
			QItemSelectionModel::Select | QItemSelectionModel::Rows);
	}
	if(not currentSelectionBackup.empty())
	{
		QModelIndex current = mdl->pathToIndex(currentSelectionBackup);
		// if the path for "current" does not exist, look for the first existing ascendant.
		while(not current.isValid() and not currentSelectionBackup.empty())
		{
			currentSelectionBackup.removeTail();	// try the parent.
			current = mdl->pathToIndex(currentSelectionBackup);
		}
		// if no ascendant exist anymore, use the configuration root.
		if(not current.isValid())
		{
			current = getConfigRoot();
		}
		selection->setCurrentIndex(current, QItemSelectionModel::Current | QItemSelectionModel::Rows);
	}
}

/**
 * Searches for an attribute matching given criteria.
 * @param what String to search for.
 * @param caseSensitive Whether the search should be case sensitive or not.
 * @param exactMatch Whether the match should simply contain the pattern or be an exact match.
 * @param column Index of the column to search in.
 */
void ConfigEditorTree::findAttribute(const QString & what, bool caseSensitive, bool exactMatch, int column)
{
	if (what.isEmpty()) return;

	// Search all occurrences, starting from configuration root by default
	QModelIndex startIndex =  model()->index(0, column, QModelIndex());
	Qt::MatchFlags flags = Qt::MatchRecursive;
	if (caseSensitive) flags |= Qt::MatchCaseSensitive;
	if (not exactMatch) flags |= Qt::MatchContains;
	QModelIndexList matches = model()->match(startIndex, Qt::DisplayRole, what, -1, flags);

	// Select next occurrence
	if (not matches.isEmpty())
	{
		qSort(matches);
		QItemSelectionModel* selection = selectionModel();
		QModelIndex currentIndex = selection->currentIndex();

		QModelIndex nextIndex = matches.first();
		if (currentIndex.isValid() and currentIndex < matches.last())
		{
			for (int i=0; i < matches.size(); ++i)
			{
				if (currentIndex < nextIndex)
					break;
				else
					nextIndex = matches[i];
			}
		}

		scrollTo(nextIndex);
		selectionModel()->select(nextIndex, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Current | QItemSelectionModel::Rows);
		selectionModel()->setCurrentIndex(nextIndex, QItemSelectionModel::Current | QItemSelectionModel::Rows);
	}
}

}
}
