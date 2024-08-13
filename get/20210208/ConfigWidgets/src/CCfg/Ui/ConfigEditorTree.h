#ifndef CCFG_UI_CONFIGEDITORTREE_H_INCLUDED
#define CCFG_UI_CONFIGEDITORTREE_H_INCLUDED

#include "../Model/AbstractConfigModelAccessor.h" // for typedef IndexList and ConfigNode
#include "CCfg/Path.h"
#include <QTreeView>
#include <set>
#include <list>
#include <string>
#include <utility>

namespace CCfg
{
using Model::ConfigNode;
using Model::IndexList;

namespace Ui
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * A tree view enabling advanced display and edition of a configuration.
 *
 * @author fchateau
 */
class ConfigEditorTree : public QTreeView
{
	Q_OBJECT

	/// backups the path of expanded items when an update requires model
	/// recomputation, incurring rows removal and reinsertion.
	std::set<std::string> expandState;

	std::list<Path> selectionBackup;
	Path currentSelectionBackup;

public:
	ConfigEditorTree(QWidget* parent = 0);
	virtual ~ConfigEditorTree();

	virtual void setModel(QAbstractItemModel* model);

	QModelIndex getConfigRoot() const;

	void expandAll(const QModelIndex& item);
	void collapseAll(const QModelIndex& item);
	void findAttribute(const QString & searchText, bool caseSensitive=false, bool exactMatch=false, int column=0);

protected:
	virtual void rowsAboutToBeRemoved(const QModelIndex& parent, int start, int end);
	virtual void rowsInserted(const QModelIndex& parent, int start, int end);

protected slots:
	virtual void onEndInsertion(ConfigNode* parent, const IndexList& indexes);
	virtual void onBeginRemoval(ConfigNode* parent, const IndexList& indexes);

	virtual void onBeginUpdate();
	virtual void onEndUpdate();

private:
	void saveExpandingState(const QModelIndex& parentIdx);
	void restoreExpandingState(const QModelIndex& parentIdx);

	void saveSelection(const QModelIndex& parentIdx);
	void restoreSelection();
};

}
}

#endif // CCFG_UI_CONFIGEDITORTREE_H_INCLUDED
