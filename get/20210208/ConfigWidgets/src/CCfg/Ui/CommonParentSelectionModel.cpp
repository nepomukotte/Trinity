#include "CommonParentSelectionModel.h"

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::CommonParentSelectionModel class methods body
////////////////////////////////////////////////////////////////////////////////


/**
 * Constructor.
 * @param model the model on which selection operates.
 * @param parent the owner of this new object.
 */
CommonParentSelectionModel::CommonParentSelectionModel(QAbstractItemModel* model,
	QObject* parent)
	: QItemSelectionModel(model, parent)
{
}

/**
 * Destructor.
 */
CommonParentSelectionModel::~CommonParentSelectionModel()
{
}

/**
 * Selects the model item index using the specified command, and emits selectionChanged().
 */
void CommonParentSelectionModel::select(const QModelIndex& index,
	QItemSelectionModel::SelectionFlags command)
{
	QItemSelection oldSelection = selection();
	// if the selection command does not clear the old selection
	// when there is a old selection.
	if((command bitand QItemSelectionModel::Clear) == 0 and oldSelection.size() > 0)
	{
		// if the item to select does not have the same parent as items of the old
		// selection.
		if(index.parent() != oldSelection[0].parent())
		{
			return;	// => discard selection command
		}
	}
	QItemSelectionModel::select(index, command);
}

/**
 * Selects the item selection using the specified command, and emits selectionChanged().
 */
void CommonParentSelectionModel::select(const QItemSelection& newSelection,
	QItemSelectionModel::SelectionFlags command)
{
	QItemSelection filteredSelection;
	if(newSelection.size() > 0)
	{
		QItemSelection oldSelection = selection();
		QModelIndex parent;
		// if the selection command does not clear the old selection
		// when there is a old selection.
		if((command bitand QItemSelectionModel::Clear) == 0 and oldSelection.size() > 0)
		{
			// all items in the new selection will have to have the same parent
			// as the items of the old selection.
			parent = oldSelection[0].parent();
		}
		else
		{
			// all items in the new selection will have to have the same parent
			// as the first item in that new selection.
			parent =  newSelection[0].parent();
		}
		// now, we filter range where the parent is not the same as the chosen one.
		foreach(QItemSelectionRange range, newSelection)
		{
			if(range.parent() == parent)
			{
				filteredSelection.append(range);
			}
		}
	}
	// forward command on the filtered selection to the base class.
	QItemSelectionModel::select(filteredSelection, command);
}

}
}
