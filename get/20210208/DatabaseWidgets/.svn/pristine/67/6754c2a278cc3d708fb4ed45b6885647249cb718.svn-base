#include "DatabaseBrowserTree.h"
#include <QHeaderView>
#include "ConnectionPropertiesDialog.h"
#include "SchemaPropertiesDialog.h"
#include "UserAccessDialog.h"
#include "../Model/DatabaseModel.h"
#include "../Model/ConnectionItem.h"
#include "../Model/ConfigItem.h"
#include "../Model/VersionItem.h"
#include "../Model/SiteItem.h"
#include "../Model/RunItem.h"
#include "../Model/ConditionItem.h"
#include "../Model/AbstractDatabaseItem.h"
#include "Utils/Error/ExcepHandling.h"
#include <QItemSelectionModel>
#include <QInputDialog>
#include <QMessageBox>
#include <QPushButton>
#include <QMenu>
#include <QAction>

namespace CCfg
{
using namespace Model;
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::DatabaseBrowserTree class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
DatabaseBrowserTree::DatabaseBrowserTree(QWidget* parent) : QTreeView(parent), acceptMode(), filterMode(DatabaseFilterProxyModel::NoFilter)
{
	tree = this;

	header()->hide();

	setAcceptMode(AcceptLoad);

	initActions();
	initMenus();

	connect(this, SIGNAL(doubleClicked(QModelIndex)),
			this, SLOT(onDoubleClicked(QModelIndex)));
}

/**
 * Destructor.
 */
DatabaseBrowserTree::~DatabaseBrowserTree()
{
}

/**
 * Defines the model for the database browser.
 */
void DatabaseBrowserTree::setDatabaseModel(Model::DatabaseModel* model)
{
	DatabaseFilterProxyModel* filteredModel = new DatabaseFilterProxyModel(this, filterMode);
	filteredModel->setSourceModel(model);

	setModel(filteredModel);
	QItemSelectionModel* selecModel = new SelectionModel(filteredModel, *this);
	setSelectionModel(selecModel);
	connect(selecModel, SIGNAL(currentChanged(QModelIndex, QModelIndex)),
			this, SLOT(onCurrentItemChanged(QModelIndex, QModelIndex)));
	connect(selecModel, SIGNAL(selectionChanged(QItemSelection, QItemSelection)),
	        this, SLOT(onSelectionChanged(QItemSelection, QItemSelection)));
	//expand(filteredModel->index(0, 0, QModelIndex())); // expand the connections root
}

/**
 * Returns the DatabaseFilterProxyModel of the tree widget.
 */
Model::DatabaseFilterProxyModel* DatabaseBrowserTree::getFilterModel() const
{
	return static_cast<Model::DatabaseFilterProxyModel*>(model());
}

/**
 * Returns the DatabaseModel of the tree widget.
 */
Model::DatabaseModel* DatabaseBrowserTree::getDatabaseModel() const
{
	DatabaseFilterProxyModel* filterModel = getFilterModel();
	if (not filterModel) return 0;
	return static_cast<Model::DatabaseModel*>(filterModel->sourceModel());
}

/**
 * Defines whether the dialog box allows to select versions to load or a
 * configuration to save.
 */
void DatabaseBrowserTree::setAcceptMode(AcceptMode newMode)
{
	acceptMode = newMode;
	setSelectionMode(acceptMode == AcceptLoad ?
			QAbstractItemView::ExtendedSelection :
			QAbstractItemView::SingleSelection);
}

/**
 * Sets schema display mode: display config versions, run conditions or both.
 */
void DatabaseBrowserTree::setFilterMode(DatabaseFilterProxyModel::FilterMode newMode)
{
	filterMode = newMode;
	DatabaseFilterProxyModel* filterModel = getFilterModel();
	if (filterModel)
		filterModel->setFilterMode(newMode);
}

/**
 * Returns the list of selected items.
 */
std::vector<Model::AbstractDatabaseItem*> DatabaseBrowserTree::getSelection()
{
	std::vector<AbstractDatabaseItem*> selection;
	QModelIndexList proxyIndexes = selectionModel()->selectedIndexes();
	foreach(QModelIndex proxyIndex, proxyIndexes)
	{
		selection.push_back(DatabaseModel::indexToItem(getFilterModel()->mapToSource(proxyIndex)));
	}
	return selection;
}


/**
 * Returns the tree item having the focus.
 */
Model::AbstractDatabaseItem* DatabaseBrowserTree::getCurrentItem()
{
	QModelIndex proxyIndex = selectionModel()->currentIndex();
	QModelIndex index = getFilterModel()->mapToSource(proxyIndex);
	return DatabaseModel::indexToItem(index);
}

/**
 * Returns the index of the tree item having the focus.
 */
QModelIndex DatabaseBrowserTree::getCurrentIndex() const
{
	QModelIndex proxyIndex = selectionModel()->currentIndex();
	return proxyIndex;
}


/**
 * Lets the user edit the properties of a connection.
 */
void DatabaseBrowserTree::editConnection(ConnectionItem* conn)
{
	ConnectionPropertiesDialog dlg(this);
	dlg.setName(conn->getName());
	dlg.setDbmsType(conn->getDbmsType());
	dlg.setUser(conn->getUser());
	dlg.setHost(conn->getHost());
	dlg.setPort(conn->getPort());
	dlg.setSID(conn->getSID());
	dlg.setSchema(conn->getSchema());
	if(dlg.exec() == QDialog::Accepted)
	{
		bool wasConnected = conn->isConnected();
		if(wasConnected)
		{
			disconnect(conn);
		}
		conn->setName(dlg.getName());
		conn->setDbmsType(dlg.getDbmsType());
		conn->setUser(dlg.getUser());
		conn->setHost(dlg.getHost());
		conn->setPort(dlg.getPort());
		conn->setSID(dlg.getSID());
		conn->setSchema(dlg.getSchema());
		if(wasConnected)
		{
			connect(conn);
		}
	}
}

/**
 *
 */
void DatabaseBrowserTree::editSchema(SchemaItem* schema)
{
	SchemaPropertiesDialog dlg(this);
	dlg.setSchemaName(schema->getName());
	dlg.setNameEditable(false);
	if(dlg.exec() == QDialog::Accepted)
	{
		schema->changeUserPassword(dlg.getPassword());
	}
}

/**
 * Changes the name of the specified item.
 */
void DatabaseBrowserTree::rename(Model::AbstractDatabaseItem* item, const std::string& newName)
{
	item->setName(newName);
}

/**
 * Removes the specified item from the tree and from its physical representation.
 */
void DatabaseBrowserTree::remove(Model::AbstractDatabaseItem* item)
{
	if(QMessageBox::warning(this, "Confirm removal",
		"The selected item will be erased permanently. Are you sure?",
		QMessageBox::Yes | QMessageBox::Cancel,
		QMessageBox::Cancel) == QMessageBox::Yes)
	{
		DatabaseModel* model = getDatabaseModel();
		QModelIndex parent = model->itemToIndex(item->getParent());
		model->removeRows(model->getRowInParent(item), 1, parent);
	}
}

/**
 * Recomputes the content of the specified item.
 */
void DatabaseBrowserTree::refresh(Model::AbstractDatabaseItem* item)
{
	DatabaseModel* model = getDatabaseModel();
	model->refresh(model->itemToIndex(item));
}

/**
 * Connects the specified connection item to its database server.
 * Afterwards, the item is filled with the list of schema and is expanded.
 */
void DatabaseBrowserTree::connect(Model::ConnectionItem* conn)
{
	if(conn->isConnected())
	{
		disconnect(conn);
	}
	conn->connect();
	if(conn->isConnected())
	{
		DatabaseModel* model = getDatabaseModel();
		QModelIndex connIndex = model->itemToIndex(conn);
		model->refresh(connIndex);

		expand(connIndex);
	}
}

/**
 * Disconnects the specified connection item from its database server.
 * The content of the item is cleared and the item is then collapsed.
 */
void DatabaseBrowserTree::disconnect(Model::ConnectionItem* conn)
{
	DatabaseModel* model = getDatabaseModel();
	QModelIndex connIndex = model->itemToIndex(conn);
	model->clear(connIndex);
	conn->disconnect();

	collapse(connIndex);
}

/**
 * This slots is called whenever the selected database item changes.
 */
void DatabaseBrowserTree::onCurrentItemChanged(const QModelIndex& currentIndex,
		const QModelIndex& /*previous*/)
{
	try
	{
		if (not currentIndex.isValid())
		{
			emit currentDatabaseItemChanged(0);
		}
		else
		{
			AbstractDatabaseItem* item = DatabaseModel::indexToItem(getFilterModel()->mapToSource(currentIndex));
			emit currentDatabaseItemChanged(item);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a context menu request.
 */
void DatabaseBrowserTree::onContextMenu(const QPoint& pos)
{
	try
	{
		AbstractDatabaseItem* item = getCurrentItem();
		enableActions(item);
		contextMenu[item->getType()]->exec(viewport()->mapToGlobal(pos));
	}
	CATCH_AND_REPORT
}

/**
 * Accepts the dialog box when the user double clicks on configs or versions.
 */
void DatabaseBrowserTree::onDoubleClicked(const QModelIndex& index)
{
	try
	{
		AbstractDatabaseItem* item = DatabaseModel::indexToItem(getFilterModel()->mapToSource(index));
		if ((acceptMode == AcceptLoad and (item->getType() == item->Version or item->getType() == item->RunCondition)) or
			(acceptMode == AcceptSave and item->getType() == item->Config))
		{
			emit databaseItemDoubleClicked(item);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Emits signal when selected items changed.
 */
void DatabaseBrowserTree::onSelectionChanged(const QItemSelection &, const QItemSelection &)
{
	emit selectionChanged();
}

/**
 *
 */
void DatabaseBrowserTree::onEdit()
{
	try
	{
		AbstractDatabaseItem* item = getCurrentItem();
		switch(item->getType())
		{
		case AbstractDatabaseItem::Connection:
			editConnection(dynamic_cast<ConnectionItem*>(item));
			break;
		case AbstractDatabaseItem::Schema:
			editSchema(dynamic_cast<SchemaItem*>(item));
			break;
		default:
			break;
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the "Rename" action.
 */
void DatabaseBrowserTree::onRename()
{
	try
	{
		tree->edit(getCurrentIndex());
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the "Delete" action.
 */
void DatabaseBrowserTree::onDelete()
{
	try
	{
		AbstractDatabaseItem* item = getCurrentItem();
		remove(item);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the "Refresh" action.
 */
void DatabaseBrowserTree::onRefresh()
{
	try
	{
		refresh(getCurrentItem());
	}
	CATCH_AND_REPORT
}

/**
 *
 */
void DatabaseBrowserTree::onConnecCreate()
{
	try
	{
		ConnectionPropertiesDialog dlg(this);
		if(dlg.exec() == QDialog::Accepted)
		{
			getDatabaseModel()->createConnection(dlg.getName(),
					dlg.getDbmsType(),
					dlg.getHost(),
					dlg.getPort(),
					dlg.getUser(),
					dlg.getSID(),
					dlg.getSchema());
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the "Connect" action.
 */
void DatabaseBrowserTree::onConnect()
{
	try
	{
		connect(dynamic_cast<ConnectionItem*>(getCurrentItem()));
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the "Connect As" action.
 */
void DatabaseBrowserTree::onConnectAs()
{
	try
	{
		Model::ConnectionItem* conn =
			dynamic_cast<Model::ConnectionItem*>(getCurrentItem());
		conn->setUser(std::string());
		conn->setPassword(std::string());
		connect(conn);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the "Disconnect" action.
 */
void DatabaseBrowserTree::onDisconnect()
{
	try
	{
		disconnect(dynamic_cast<ConnectionItem*>(getCurrentItem()));
	}
	CATCH_AND_REPORT
}


/**
 *
 */
void DatabaseBrowserTree::onUserCreate()
{
	try
	{
		SchemaPropertiesDialog dlg(this);
		if(dlg.exec() == QDialog::Accepted)
		{
			ConnectionItem* parent = dynamic_cast<ConnectionItem*>(getCurrentItem());
			getDatabaseModel()->createSchema(parent, dlg.getSchemaName(), dlg.getPassword());
		}
	}
	CATCH_AND_REPORT
}

/**
 *
 */
void DatabaseBrowserTree::onUserGrant()
{
	try
	{
		SchemaItem* schema = dynamic_cast<SchemaItem*>(getCurrentItem());
		UserAccessDialog dlg(schema->getDao(), this);
		dlg.setSchema(schema->getName());
		dlg.exec();
	}
	CATCH_AND_REPORT
}

/**
 *
 */
void DatabaseBrowserTree::onConfigCreate()
{
	try
	{
		SchemaItem* schema = dynamic_cast<SchemaItem*>(getCurrentItem());
		bool ok;
		QString name = QInputDialog::getText(this, tr("New Configuration"),
				tr("Config name:"), QLineEdit::Normal, QString(), &ok);
		if(ok and not name.isEmpty())
		{
			getDatabaseModel()->createConfig(schema, name.toStdString());
		}
	}
	CATCH_AND_REPORT
}

/**
 *	Creates a new site in the database for storing run conditions.
 */
void DatabaseBrowserTree::onSiteCreate()
{
	try
	{
		SchemaItem* schema = dynamic_cast<SchemaItem*>(getCurrentItem());
		bool ok;
		QString name = QInputDialog::getText(this, tr("New Site"),
				tr("Site name:"), QLineEdit::Normal, QString(), &ok);
		if(ok and not name.isEmpty())
		{
			getDatabaseModel()->createSite(schema, name.toStdString());
		}
	}
	CATCH_AND_REPORT
}

/**
 *	Creates a new run for a site in the database.
 */
void DatabaseBrowserTree::onRunCreate()
{
	try
	{
		SiteItem* site = dynamic_cast<SiteItem*>(getCurrentItem());
		getDatabaseModel()->createRun(site);
	}
	CATCH_AND_REPORT
}


/**
 * Initializes the actions of the contextual menues.
 */
void DatabaseBrowserTree::initActions()
{
	actionEdit 		= new QAction(tr("Edit..."), this);
	actionRename	= new QAction(tr("Rename"), this);
	actionDelete	= new QAction(tr("Delete"), this);
	actionRefresh	= new QAction(tr("Refresh"), this);

	actionConnecCreate	= new QAction(tr("Create connection..."), this);
	actionConnect		= new QAction(tr("Connect"), this);
	actionConnectAs		= new QAction(tr("Connect as..."), this);
	actionDisconnect	= new QAction(tr("Disconnect"), this);
	actionUserCreate	= new QAction(tr("Create user/schema..."), this);
	actionUserGrant		= new QAction(tr("Edit access..."), this);
	actionConfigCreate	= new QAction(tr("Create config..."), this);
	actionSiteCreate	= new QAction(tr("Create site..."), this);
	actionRunCreate		= new QAction(tr("Create run..."), this);

	connect(actionEdit,			SIGNAL(triggered()), this, SLOT(onEdit()));
	connect(actionRename,		SIGNAL(triggered()), this, SLOT(onRename()));
	connect(actionDelete,		SIGNAL(triggered()), this, SLOT(onDelete()));
	connect(actionRefresh,		SIGNAL(triggered()), this, SLOT(onRefresh()));
	connect(actionConnecCreate, SIGNAL(triggered()), this, SLOT(onConnecCreate()));
	connect(actionConnect,		SIGNAL(triggered()), this, SLOT(onConnect()));
	connect(actionConnectAs,	SIGNAL(triggered()), this, SLOT(onConnectAs()));
	connect(actionDisconnect,	SIGNAL(triggered()), this, SLOT(onDisconnect()));
	connect(actionUserCreate,	SIGNAL(triggered()), this, SLOT(onUserCreate()));
	connect(actionUserGrant,	SIGNAL(triggered()), this, SLOT(onUserGrant()));
	connect(actionConfigCreate, SIGNAL(triggered()), this, SLOT(onConfigCreate()));
	connect(actionSiteCreate, 	SIGNAL(triggered()), this, SLOT(onSiteCreate()));
	connect(actionRunCreate, 	SIGNAL(triggered()), this, SLOT(onRunCreate()));
}

/**
 * Initializes the contextual menues.
 */
void DatabaseBrowserTree::initMenus()
{
	QMenu*& menuConnecRoot	= contextMenu[AbstractDatabaseItem::Root];
	QMenu*& menuConnec		= contextMenu[AbstractDatabaseItem::Connection];
	QMenu*& menuUser		= contextMenu[AbstractDatabaseItem::Schema];
	QMenu*& menuConfig		= contextMenu[AbstractDatabaseItem::Config];
	QMenu*& menuVersion		= contextMenu[AbstractDatabaseItem::Version];
	QMenu*& menuSite		= contextMenu[AbstractDatabaseItem::Site];
	QMenu*& menuRun			= contextMenu[AbstractDatabaseItem::Run];
	QMenu*& menuCondition	= contextMenu[AbstractDatabaseItem::RunCondition];

	menuConnecRoot	= new QMenu(this);
	menuConnec		= new QMenu(this);
	menuUser		= new QMenu(this);
	menuConfig		= new QMenu(this);
	menuVersion		= new QMenu(this);
	menuSite		= new QMenu(this);
	menuRun			= new QMenu(this);
	menuCondition	= new QMenu(this);

	menuConnecRoot->addAction(actionConnecCreate);

	menuConnec->addAction(actionEdit);
	menuConnec->addAction(actionRename);
	menuConnec->addAction(actionDelete);
	menuConnec->addAction(actionRefresh);
	menuConnec->addSeparator();
	menuConnec->addAction(actionConnect);
	menuConnec->addAction(actionConnectAs);
	menuConnec->addAction(actionDisconnect);
	menuConnec->addSeparator();
	menuConnec->addAction(actionUserCreate);

	menuUser->addAction(actionEdit);
	menuUser->addAction(actionDelete);
	menuUser->addAction(actionRefresh);
	menuUser->addSeparator();
	menuUser->addAction(actionUserGrant);
	menuUser->addAction(actionConfigCreate);
	menuUser->addAction(actionSiteCreate);

	menuConfig->addAction(actionRename);
	menuConfig->addAction(actionDelete);
	menuConfig->addAction(actionRefresh);

	menuVersion->addAction(actionRename);
	menuVersion->addAction(actionDelete);

	menuSite->addAction(actionRename);
	menuSite->addAction(actionRunCreate);
	menuSite->addAction(actionDelete);

	menuRun->addAction(actionDelete);

	menuCondition->addAction(actionDelete);

	tree->setContextMenuPolicy(Qt::CustomContextMenu);
	connect(tree, SIGNAL(customContextMenuRequested(QPoint)),
		this, SLOT(onContextMenu(QPoint)));
}

/**
 * Enables or disable the context menu actions according to the state of the item.
 */
void DatabaseBrowserTree::enableActions(Model::AbstractDatabaseItem* item)
{
	ConnectionItem* conn = dynamic_cast<ConnectionItem*>(item);
	if(conn != 0)
	{
		bool isConnected = conn->isConnected();
		actionUserCreate->setEnabled(isConnected);
		actionDisconnect->setEnabled(isConnected);
	}
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::DatabaseBrowserTree::SelectionModel class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
DatabaseBrowserTree::SelectionModel::SelectionModel(QAbstractItemModel* model, DatabaseBrowserTree& browser)
	: QItemSelectionModel(model, &browser), browser(browser)
{
}

/**
 * Destructor.
 */
DatabaseBrowserTree::SelectionModel::~SelectionModel()
{
}

/*
 * (see inherited doc)
 */
void DatabaseBrowserTree::SelectionModel::select(const QModelIndex& index, QItemSelectionModel::SelectionFlags command)
{
	if(isSelectable(index))
	{
		QItemSelectionModel::select(index, command);
	}
}

/*
 * (see inherited doc)
 */
void DatabaseBrowserTree::SelectionModel::select(const QItemSelection& selection, QItemSelectionModel::SelectionFlags command)
{
	QItemSelection filteredSelection;
	if(selection.size() > 0)
	{
		foreach(QItemSelectionRange range, selection)
		{
			if(isSelectable(range.topLeft()))
			{
				filteredSelection.append(range);
			}
		}
	}
	QItemSelectionModel::select(filteredSelection, command);
}

/**
 * Returns whether the specified item can be selected in the current dialog mode.
 */
bool DatabaseBrowserTree::SelectionModel::isSelectable(const QModelIndex& proxyIndex) const
{
	QModelIndex index = browser.getFilterModel()->mapToSource(proxyIndex);
	if(not index.isValid())	return false;
	AbstractDatabaseItem* item = DatabaseModel::indexToItem(index);
	AbstractDatabaseItem::ItemType type = item->getType();
	return
		(browser.getAcceptMode() == AcceptLoad
			and (type == AbstractDatabaseItem::Version or type == AbstractDatabaseItem::RunCondition))
		or
		(browser.getAcceptMode() == AcceptSave
			and (type == AbstractDatabaseItem::Config or type == AbstractDatabaseItem::Run));
}

}
}
