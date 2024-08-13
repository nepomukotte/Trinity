#ifndef CCFG_UI_DATABASEBROWSERTREE_H_INCLUDED
#define CCFG_UI_DATABASEBROWSERTREE_H_INCLUDED

#include <QTreeView>
#include "CCfg/Model/AbstractDatabaseItem.h"
#include "../Model/SchemaItem.h"
#include "../Model/DatabaseFilterProxyModel.h"
#include <QItemSelectionModel>
#include <QModelIndex>
#include <vector>
#include <memory>

// forward declarations
class QModelIndex;
class QAction;
class QMenu;

namespace CCfg
{
// forward declarations
namespace Model
{
class DatabaseModel;
class ConnectionItem;
class ConfigItem;
class VersionItem;
}
using Model::AbstractDatabaseItem;
using Model::SchemaItem;
using Model::DatabaseFilterProxyModel;
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 *
 * @author fchateau
 */
class DatabaseBrowserTree : public QTreeView
{
	Q_OBJECT
public:
	// enums
	enum AcceptMode { AcceptLoad, AcceptSave };

	// classes
	class SelectionModel;

	// methods
	DatabaseBrowserTree(QWidget* parent = 0);
	virtual ~DatabaseBrowserTree();

	void setDatabaseModel(Model::DatabaseModel* model);
	void setAcceptMode(AcceptMode newMode);
	Model::DatabaseModel* getDatabaseModel() const;
	AcceptMode getAcceptMode() const { return acceptMode; }
	Model::DatabaseFilterProxyModel* getFilterModel() const;
	void setFilterMode(DatabaseFilterProxyModel::FilterMode newMode);

	std::vector<Model::AbstractDatabaseItem*> getSelection();
	Model::AbstractDatabaseItem* getCurrentItem();
	QModelIndex getCurrentIndex() const;
signals:
	void currentDatabaseItemChanged(AbstractDatabaseItem* newCurrrentItem);
	void databaseItemDoubleClicked(AbstractDatabaseItem* item);
	void selectionChanged();

protected:
	void editConnection(Model::ConnectionItem* conn);
	void editSchema(Model::SchemaItem* schema);
	void rename(Model::AbstractDatabaseItem* item, const std::string& newName);
	void remove(Model::AbstractDatabaseItem* item);
	void refresh(Model::AbstractDatabaseItem* item);
	void connect(Model::ConnectionItem* conn);
	void disconnect(Model::ConnectionItem* conn);

	using QObject::connect;

protected slots:
	void onCurrentItemChanged(const QModelIndex& current, const QModelIndex& previous);
	void onContextMenu(const QPoint& pos);
	void onDoubleClicked(const QModelIndex& index);
	void onSelectionChanged(const QItemSelection &, const QItemSelection &);

	void onEdit();
	void onRename();
	void onDelete();
	void onRefresh();
	void onConnecCreate();
	void onConnect();
	void onConnectAs();
	void onDisconnect();
	void onUserCreate();
	void onUserGrant();
	void onConfigCreate();
	void onSiteCreate();
	void onRunCreate();

private:
	void initActions();
	void initMenus();
	void enableActions(Model::AbstractDatabaseItem* item);

	DatabaseBrowserTree(const DatabaseBrowserTree& rvalue);
	DatabaseBrowserTree& operator=(const DatabaseBrowserTree& rvalue);

	AcceptMode acceptMode;
	DatabaseFilterProxyModel::FilterMode filterMode;

	QMenu* contextMenu[Model::AbstractDatabaseItem::TypeCount];

	QAction* actionEdit;
	QAction* actionRename;
	QAction* actionDelete;
	QAction* actionRefresh;

	QAction* actionConnecCreate;
	QAction* actionConnect;
	QAction* actionConnectAs;
	QAction* actionDisconnect;
	QAction* actionUserCreate;
	QAction* actionUserGrant;
	QAction* actionConfigCreate;
	QAction* actionSiteCreate;
	QAction* actionRunCreate;

	QTreeView* tree;
};

/**
 * This class enforces selection restrictions on the database browser tree.
 * Only version items can be selected when the dialog is in AcceptLoad mode and
 * only config items when in AcceptSave mode.
 */
class DatabaseBrowserTree::SelectionModel : public QItemSelectionModel
{
	Q_OBJECT
	DatabaseBrowserTree& browser;
public:
	SelectionModel(QAbstractItemModel* model, DatabaseBrowserTree& browser);
	~SelectionModel();

public slots:
	void select(const QModelIndex& index, QItemSelectionModel::SelectionFlags command);
	void select(const QItemSelection& selection, QItemSelectionModel::SelectionFlags command);

protected:
	bool isSelectable(const QModelIndex& index) const;
};

}
}

#endif // CCFG_UI_DATABASEBROWSERTREE_H_INCLUDED
