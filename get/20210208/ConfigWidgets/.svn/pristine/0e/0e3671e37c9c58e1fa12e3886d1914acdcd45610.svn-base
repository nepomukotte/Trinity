#include "ConfigEditorController.h"
#include "Delegate/ConfigEditorDelegate.h"
#include "../Ui/ConfigEditorTree.h"
#include "../Model/ConfigModel.h"
#include "../Model/ConfigNode.h"
#include "../Model/AttributeBuilderModel.h"
#include "Utils/Error/ExcepHandling.h"
#include "CCfg/Io/Url.h"
#include "CCfg/Attribute.h"
#include <QMenu>
#include <QDialog>
#include <QPoint>
#include <QModelIndex>
#include <QItemSelection>
#include <QMessageBox>
#include <QInputDialog>
#include <QFileDialog>

namespace CCfg
{
using Model::ConfigModel;
using Model::ConfigNode;

namespace Control
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Control::ConfigEditorController class methods body
///////////////////////////////////////////////////////////////////////////////

const char ConfigEditorController::FILE_EXTENSIONS_FILTER[] =
	"XML config files (*.xml *.xcfg)";

/**
 * Default constructor.
 * @param parent the parent QObject.
 */
ConfigEditorController::ConfigEditorController(QObject* parent)
	: BaseConfigEditorController(parent)
{
}

/**
 * Destructor.
 */
ConfigEditorController::~ConfigEditorController()
{
}

/*
 * doc: see base overriden method
 */
void ConfigEditorController::setView(Ui::ConfigEditorTree* newView)
{
	if(newView == view)
	{
		return;
	}
	if(view != 0)
	{
		disconnect(view, 0,	this, 0);
	}
	BaseConfigEditorController::setView(newView);

	if(newView != 0)
	{
		Delegate::ConfigEditorDelegate* delegate = new Delegate::ConfigEditorDelegate(newView);
		delegate->setFileDialogFilter(ConfigEditorController::FILE_EXTENSIONS_FILTER);
		setDelegate(delegate);

		initActions(newView);
		initMenus();

		newView->setContextMenuPolicy(Qt::CustomContextMenu);
		connect(newView, SIGNAL(customContextMenuRequested(QPoint)),
			this, SLOT(onContextMenu(QPoint)));
	}
}

/*
 *
 */
void ConfigEditorController::onInsertNew(Model::AttributeBuilderModel* builder,
	const QString& /*path*/, bool insertAsChild)
{
	try
	{
		std::auto_ptr<Attribute> newAttribute(builder->createAttribute());
		QItemSelectionModel* selection = view->selectionModel();

		if(insertAsChild)
		{
			QModelIndex parent = selection->currentIndex();
			model->insertInside(parent, newAttribute);
		}
		else
		{
			QModelIndex sibling = selection->currentIndex();
			model->insertAfter(sibling, newAttribute);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the user request of a context menu on the view.
 * @param pos the point (in the view's viewport coordinates) where the user
 * invoked the menu.
 */
void ConfigEditorController::onContextMenu(const QPoint& pos)
{
	try
	{
		QMenu *menu = 0;
		if(isMultiselection(view->selectionModel()->selection()))
		{
			menu = menuMultiselection.get();
		}
		else
		{
			QModelIndexList indexes = view->selectionModel()->selectedIndexes();
			assert(indexes.size() <= 1);
			if(indexes.size() == 1)
			{
				ConfigNode* node = model->indexToNode(indexes.front());
				if(node->isRoot())
				{
					menu = menuRootNode.get();
				}
				else
				{
					if(node->isObject())
					{
						menu = menuObjectNode.get();
					}
					else
					{
						menu = menuScalarNode.get();
					}
				}
			}
		}
		if(menu != 0)
		{
			enableMenuItems();
			menu->exec(view->viewport()->mapToGlobal(pos));
		}
	}
	CATCH_AND_REPORT
}

/**
 * Computes the content of each inheritance menu when they parent menu is about to show.
 */
void ConfigEditorController::onInheritanceMenuAboutToShow()
{
	try
	{
		QModelIndexList indexes = view->selectionModel()->selectedIndexes();
		if(indexes.size() == 1)
		{
			ConfigNode* node = model->indexToNode(indexes.front());

			fillMenuOverrides(node);
			fillMenuOverridenBy(node);
			fillMenuInherits(node);
			fillMenuInheritedBy(node);
			fillMenuDefinedIn(node);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a click on a ConfigNode reference menu item.
 * It puts the selection on the node specified in the menu.
 */
void ConfigEditorController::onNodeRefClicked()
{
	try
	{
		QAction* action = qobject_cast<QAction*>(sender());
		QModelIndex index = model->pathToIndex(Path(action->text().toStdString()));
		view->selectionModel()->select(index, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
		view->selectionModel()->setCurrentIndex(index, QItemSelectionModel::Current | QItemSelectionModel::Rows);
	}
	CATCH_AND_REPORT
}

/**
 * Displays the "New attribute" dialog box.
 */
void ConfigEditorController::onNew()
{
	try
	{
		if(not isActionAuthorized(ActionNew, view->selectionModel()->selection()))
		{
			notifyUnauthorizedAction();
		}
		else
		{
			QPoint pt;
			QRect rect = view->visualRect(view->selectionModel()->currentIndex());
			if(rect.isValid())
			{
				pt = view->viewport()->mapToGlobal(rect.bottomRight());
			}
			else
			{
				pt = view->parentWidget()->mapToGlobal(view->geometry().center());
			}

			emit attributeBuilderRequested(pt);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Overrides an inherited attribute (create a local copy...).
 */
void ConfigEditorController::onOverride()
{
	try
	{
		if(not isActionAuthorized(ActionOverride, view->selectionModel()->selection()))
		{
			notifyUnauthorizedAction();
		}
		else
		{
			model->override(view->selectionModel()->currentIndex());
		}
	}
	CATCH_AND_REPORT
}

/**
 * Displays the "rename" dialog box.
 */
void ConfigEditorController::onRename()
{
	try
	{
		if(not isActionAuthorized(ActionRename, view->selectionModel()->selection()))
		{
			notifyUnauthorizedAction();
		}
		else
		{
			QModelIndex currentIndex = view->selectionModel()->currentIndex();
			ConfigNode* currentNode = model->indexToNode(currentIndex);

			bool ok;
			QString newName = QInputDialog::getText(view, tr("Rename dialog"),
					tr("New name:"), QLineEdit::Normal,
					QString::fromStdString(currentNode->getName()), &ok);
			if(ok)
			{
				model->setName(currentIndex, newName.toStdString());
			}
		}
	}
	CATCH_AND_REPORT
}

/**
 * Deletes the selected attributes from the model.
 */
void ConfigEditorController::onDelete()
{
	try
	{
		if(not isActionAuthorized(ActionDelete, view->selectionModel()->selection()))
		{
			notifyUnauthorizedAction();
		}
		else
		{
			QModelIndexList modelIndexes = view->selectionModel()->selectedRows();
			if(modelIndexes.size() > 0)
			{
				IndexList indexes;
				foreach(QModelIndex modelIndex, modelIndexes)
				{
					indexes.push_back(modelIndex.row());
				}
				model->remove(modelIndexes[0].parent(), indexes);
			}
		}
	}
	CATCH_AND_REPORT
}

/**
 * Copies the selected attributes to the clipboard.
 */
void ConfigEditorController::onCopy()
{
	try
	{
	}
	CATCH_AND_REPORT
}

/**
 * Cut the selected attributes to the clipboard.
 */
void ConfigEditorController::onCut()
{
	try
	{
	}
	CATCH_AND_REPORT
}

/**
 * Appends attributes in the clipboard to the selected object.
 */
void ConfigEditorController::onPaste()
{
	try
	{
	}
	CATCH_AND_REPORT
}

/**
 * Displays the "search" dialog box.
 */
void ConfigEditorController::onSearch()
{
	try
	{
	}
	CATCH_AND_REPORT
}

/**
 * Recursively expands all items under the current.
 */
void ConfigEditorController::onExpand()
{
	try
	{
		if(not isActionAuthorized(ActionExpand, view->selectionModel()->selection()))
		{
			notifyUnauthorizedAction();
		}
		else
		{
			view->expandAll(view->selectionModel()->currentIndex());
		}
	}
	CATCH_AND_REPORT
}

/**
 * Recursively collapses all items under the current.
 */
void ConfigEditorController::onCollapse()
{
	try
	{
		if(not isActionAuthorized(ActionCollapse, view->selectionModel()->selection()))
		{
			notifyUnauthorizedAction();
		}
		else
		{
			view->collapseAll(view->selectionModel()->currentIndex());
		}
	}
	CATCH_AND_REPORT
}

/**
 * Appends the content of a configuration file to the selected object.
 */
void ConfigEditorController::onImportFile()
{
	try
	{
		QString filePath = selectInputFile();
		if(not filePath.isEmpty())
		{
			QModelIndex currentIndex = view->selectionModel()->currentIndex();

			QDir file(filePath);
			std::string url = Io::Url::PathNormalizer::escape(file.absolutePath().toStdString());
			model->import(currentIndex, url);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Appends the content of a configuration database entry to the selected object.
 */
void ConfigEditorController::onImportDB()
{
	try
	{
	}
	CATCH_AND_REPORT
}

/**
 * Writes the selected object and its subtree in a configuration file.
 */
void ConfigEditorController::onExportFile()
{
	try
	{
		QString filePath = selectOutputFile();
		if(not filePath.isEmpty())
		{
			QModelIndex currentIndex = view->selectionModel()->currentIndex();

			QDir file(filePath);
			std::string url = Io::Url::PathNormalizer::escape(file.absolutePath().toStdString());
			model->exportSubConfig(currentIndex, url);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Writes the selected object and its subtree in a configuration database.
 */
void ConfigEditorController::onExportDB()
{
	try
	{
	}
	CATCH_AND_REPORT
}

/**
 * Links the content of the specified configuration file to the selected object.
 */
void ConfigEditorController::onIncludeFile()
{
	try
	{
		QString filePath = selectInputFile();
		if(not filePath.isEmpty())
		{
			QModelIndex currentIndex = view->selectionModel()->currentIndex();

			QDir nodeDir = ConfigModel::getConfigDirectory(currentIndex);
			QString relativeFilePath = nodeDir.relativeFilePath(filePath);
			std::string url = Io::Url::PathNormalizer::escape(relativeFilePath.toStdString());
			model->include(currentIndex, url);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Links the content of the specified configuration database entry to the selected object.
 */
void ConfigEditorController::onIncludeDB()
{
	try
	{
	}
	CATCH_AND_REPORT
}

/**
 * Saves the selected configuration root to its file.
 */
void ConfigEditorController::onSave()
{
	try
	{
	}
	CATCH_AND_REPORT
}

/**
 * Commits the selected configuration root to its database.
 */
void ConfigEditorController::onCommit()
{
	try
	{
	}
	CATCH_AND_REPORT
}

/**
 * Initializes contextual menus actions.
 * @param parent the parent widget for actions (shortcuts are enabled if the
 * widget has the focus).
 */
void ConfigEditorController::initActions(QWidget* parent)
{
	actionNew			= new QAction(tr("New..."), parent);
	actionOverride		= new QAction(tr("Override"), parent);
	actionRename		= new QAction(tr("Rename..."), parent);
	actionDelete		= new QAction(tr("Delete"), parent);
	actionCopy			= new QAction(tr("Copy"), parent);
	actionCut			= new QAction(tr("Cut"), parent);
	actionPaste			= new QAction(tr("Paste"), parent);
	actionSearch		= new QAction(tr("Search..."), parent);
	actionExpand		= new QAction(tr("Expand all"), parent);
	actionCollapse		= new QAction(tr("Collapse all"), parent);
	actionImportFile	= new QAction(tr("File..."), parent);
	actionImportDB		= new QAction(tr("Database..."), parent);
	actionExportFile	= new QAction(tr("File..."), parent);
	actionExportDB		= new QAction(tr("Database..."), parent);
	actionIncludeFile	= new QAction(tr("File..."), parent);
	actionIncludeDB		= new QAction(tr("Database..."), parent);

	actionDelete->setShortcut(QKeySequence::Delete);
	actionDelete->setShortcutContext(Qt::WidgetShortcut);
	parent->addAction(actionDelete);

	connect(actionNew,			SIGNAL(triggered()), this, SLOT(onNew()));
	connect(actionOverride,		SIGNAL(triggered()), this, SLOT(onOverride()));
	connect(actionRename,		SIGNAL(triggered()), this, SLOT(onRename()));
	connect(actionDelete,		SIGNAL(triggered()), this, SLOT(onDelete()));
	connect(actionCopy,			SIGNAL(triggered()), this, SLOT(onCopy()));
	connect(actionCut,			SIGNAL(triggered()), this, SLOT(onCut()));
	connect(actionPaste,		SIGNAL(triggered()), this, SLOT(onPaste()));
	connect(actionSearch,		SIGNAL(triggered()), this, SLOT(onSearch()));
	connect(actionExpand,		SIGNAL(triggered()), this, SLOT(onExpand()));
	connect(actionCollapse,		SIGNAL(triggered()), this, SLOT(onCollapse()));
	connect(actionImportFile,	SIGNAL(triggered()), this, SLOT(onImportFile()));
	connect(actionImportDB,		SIGNAL(triggered()), this, SLOT(onImportDB()));
	connect(actionExportFile,	SIGNAL(triggered()), this, SLOT(onExportFile()));
	connect(actionExportDB,		SIGNAL(triggered()), this, SLOT(onExportDB()));
	connect(actionIncludeFile,	SIGNAL(triggered()), this, SLOT(onIncludeFile()));
	connect(actionIncludeDB,	SIGNAL(triggered()), this, SLOT(onIncludeDB()));
}

/**
 * Initializes contextual menus and submenus.
 */
void ConfigEditorController::initMenus()
{
	menuRootNode.reset(new QMenu);
	menuObjectNode.reset(new QMenu);
	menuScalarNode.reset(new QMenu);
	menuMultiselection.reset(new QMenu);

	menuInheritance	= new QMenu(tr("Inheritance"), menuObjectNode.get());
	connect(menuInheritance, SIGNAL(aboutToShow()), this, SLOT(onInheritanceMenuAboutToShow()));

	menuOverrides	= new QMenu(tr("Overrides"), menuInheritance);
	menuOverridenBy	= new QMenu(tr("Overriden by"), menuInheritance);
	menuInherits	= new QMenu(tr("Inherits"), menuInheritance);
	menuInheritedBy	= new QMenu(tr("Inherited by"), menuInheritance);
	menuDefinedIn	= new QMenu(tr("Defined in"), menuInheritance);

	menuInheritance->addMenu(menuOverrides);
	menuInheritance->addMenu(menuOverridenBy);
	menuInheritance->addMenu(menuInherits);
	menuInheritance->addMenu(menuInheritedBy);
	menuInheritance->addMenu(menuDefinedIn);

	menuImport		= new QMenu(tr("Import"), menuObjectNode.get());
	menuExport		= new QMenu(tr("Export"), menuObjectNode.get());
	menuInclude		= new QMenu(tr("Include"), menuObjectNode.get());

	menuImport->addAction(actionImportFile);
	menuImport->addAction(actionImportDB);
	menuExport->addAction(actionExportFile);
	menuExport->addAction(actionExportDB);
	menuInclude->addAction(actionIncludeFile);
	menuInclude->addAction(actionIncludeDB);

	createRootContextMenu(*menuRootNode);
	createMultiselectionContextMenu(*menuMultiselection);
	createAttributeContextMenu(*menuObjectNode, true);
	createAttributeContextMenu(*menuScalarNode, false);
}

/**
 * Fills the "Overrides" submenu with an action for each attribute overriden by
 * the specified node.
 */
void ConfigEditorController::fillMenuOverrides(Model::ConfigNode* node)
{
	if(not node->isInherited())
	{
		std::vector<ConfigNode*> overriden = node->getOverriden();
		fillInheritanceMenu(*menuOverrides, overriden);
	}
	else
	{
		menuOverrides->setEnabled(false);
	}
}

/**
 * Fills the "Overriden by" submenu with an action for each attribute that overrides
 * the specified node.
 */
void ConfigEditorController::fillMenuOverridenBy(Model::ConfigNode* node)
{
	if(not node->isInherited())
	{
		std::vector<ConfigNode*> overriding = node->getOverriding();
		fillInheritanceMenu(*menuOverridenBy, overriding);
	}
	else
	{
		menuOverridenBy->setEnabled(false);
	}
}

/**
 * Fills the "Inherits" submenu with an action for each attribute inherited by
 * the specified node.
 */
void ConfigEditorController::fillMenuInherits(Model::ConfigNode* node)
{
	if(not node->isInherited())
	{
		std::vector<ConfigNode*> inherited = node->getInherited();
		fillInheritanceMenu(*menuInherits, inherited);
	}
	else
	{
		menuInherits->setEnabled(false);
	}
}

/**
 * Fills the "Inherited by" submenu with an action for each attribute that inherits
 * from the specified node.
 */
void ConfigEditorController::fillMenuInheritedBy(Model::ConfigNode* node)
{
	if(not node->isInherited())
	{
		std::vector<ConfigNode*> inheriting = node->getInheriting();
		fillInheritanceMenu(*menuInheritedBy, inheriting);
	}
	else
	{
		menuInheritedBy->setEnabled(false);
	}
}

/**
 * Fills the "Defined in" submenu with an action for the inheritance definition
 * of the specified node.
 */
void ConfigEditorController::fillMenuDefinedIn(Model::ConfigNode* node)
{
	if(node->isInherited())
	{
		std::vector<ConfigNode*> definitionList;
		ConfigNode* definition = node->getDefinition();
		if(definition != 0)
		{
			definitionList.push_back(definition);
		}
		fillInheritanceMenu(*menuDefinedIn, definitionList);
	}
	else
	{
		menuDefinedIn->setEnabled(false);
	}
}

/**
 * Enables or disables menu items according to current selection.
 */
void ConfigEditorController::enableMenuItems()
{
	const QItemSelection selection = view->selectionModel()->selection();

	actionNew->setEnabled(isActionAuthorized(ActionNew, selection));
	actionOverride->setEnabled(isActionAuthorized(ActionOverride, selection));
	actionRename->setEnabled(isActionAuthorized(ActionRename, selection));
	actionDelete->setEnabled(isActionAuthorized(ActionDelete, selection));

	actionCopy->setEnabled(isActionAuthorized(ActionCopy, selection));
	actionCut->setEnabled(isActionAuthorized(ActionCut, selection));
	actionPaste->setEnabled(isActionAuthorized(ActionPaste, selection));

	menuImport->setEnabled(isActionAuthorized(ActionImport, selection));
	menuExport->setEnabled(isActionAuthorized(ActionExport, selection));
	menuInclude->setEnabled(isActionAuthorized(ActionInclude, selection));

	actionSearch->setEnabled(isActionAuthorized(ActionSearch, selection));
	actionExpand->setEnabled(isActionAuthorized(ActionExpand, selection));
	actionCollapse->setEnabled(isActionAuthorized(ActionCollapse, selection));
}

/**
 * Determines whether an action is authorized on a given selection.
 * @param actionId the id of the action.
 * @param selection the selection.
 * @return Returns true if the action can be executed on the selection.
 */
bool ConfigEditorController::isActionAuthorized(ActionId actionId, const QItemSelection& selection) const
{
	bool authorized = true;
	int selectedCount = 0;
	// for each selection range.
	for(QItemSelection::const_iterator iRange = selection.begin(); iRange != selection.end(); ++iRange)
	{
		// for each item in the range.
		for(int row = iRange->top(); row <= iRange->bottom(); ++row)
		{
			++selectedCount;
			// the action must be authorized for each selected item.
			authorized = authorized and isActionAuthorized(actionId, model->index(row, 0, iRange->parent()));
			if(selectedCount > 1)  // if multiselection
			{
				// allows only delete, copy, and cut
				if(actionId != ActionDelete and
					actionId != ActionCopy and
					actionId != ActionCut)
				{
					return false;
				}
			}
		}
	}
	return authorized;
}

/**
 * Determines whether an action is authorized on a given node.
 * @param actionId the id of the action.
 * @param nodeIndex the index of the node.
 * @return Returns true if the action can be executed on the node.
 */
bool ConfigEditorController::isActionAuthorized(ActionId actionId, const QModelIndex& nodeIndex) const
{
	ConfigNode* node = model->indexToNode(nodeIndex);
	// copy, export and rename are always authorized
	if(actionId == ActionCopy or
		actionId == ActionExport or
		actionId == ActionRename)
	{
		return true;
	}
	else
	{
		if(node->isInherited())
		{
			return actionId == ActionOverride;
		}
		else
		{
			if(node->getParent() == 0) // root node
			{
				// allow everything except delete and cut.
				return actionId != ActionDelete and actionId != ActionCut;
			}
			else
			{
				return actionId != ActionOverride;
			}
		}
	}
}

/**
 * Notifies the user that an action he was going to execute is not authorized.
 */
void ConfigEditorController::notifyUnauthorizedAction() const
{
	QMessageBox::information(view, "Action failure",
		"<b>This action is not authorized on the current selection.</b>");
}

QString ConfigEditorController::selectInputFile()
{
	return QFileDialog::getOpenFileName(view,
	 		QString(),
	 		QDir::currentPath(),
	 		FILE_EXTENSIONS_FILTER);
}

QString ConfigEditorController::selectOutputFile()
{
	return QFileDialog::getSaveFileName(view,
	 		QString(),
	 		QDir::currentPath(),
	 		ConfigEditorController::FILE_EXTENSIONS_FILTER);
}

/**
 * Builds the context menu for the root node.
 * @param menu the menu that will host menu items.
 */
void ConfigEditorController::createRootContextMenu(QMenu& menu)
{
	menu.addAction(actionNew);
	menu.addAction(actionRename);
	menu.addSeparator();
	menu.addAction(actionCopy);
	menu.addAction(actionPaste);
	menu.addSeparator();
	menu.addAction(actionSearch);
	menu.addAction(actionExpand);
	menu.addAction(actionCollapse);
	menu.addSeparator();
	menu.addMenu(menuImport);
	menu.addMenu(menuExport);
}

/**
 * Builds the context menu for multiselection.
 * @param menu the menu that will host menu items.
 */
void ConfigEditorController::createMultiselectionContextMenu(QMenu& menu)
{
	menu.addAction(actionCopy);
	menu.addAction(actionCut);
	menu.addSeparator();
	menu.addAction(actionDelete);
}

/**
 * Builds context menus for attribute nodes.
 * @param menu the menu that will host menu items.
 * @param isObject indicates whether we request to build the contextual menu
 * for objects or for scalars.
 */
void ConfigEditorController::createAttributeContextMenu(QMenu& menu, bool isObject)
{
	menu.addAction(actionNew);
	menu.addAction(actionOverride);
	menu.addAction(actionRename);
	menu.addAction(actionDelete);
	menu.addSeparator();
	menu.addAction(actionCopy);
	menu.addAction(actionCut);
	menu.addAction(actionPaste);

	if(isObject)
	{
		menu.addSeparator();
		menu.addMenu(menuImport);
		menu.addMenu(menuExport);
		menu.addMenu(menuInclude);
		menu.addSeparator();
		menu.addAction(actionSearch);
		menu.addAction(actionExpand);
		menu.addAction(actionCollapse);
	}

	menu.addSeparator();
	menu.addMenu(menuInheritance);
}

/**
 * Fills either of the "Inherited" or "Inheriting" submenus using the specified
 * inheritance list.
 * @param menu the menu that will host created actions.
 * @param list the list of configuration nodes that must be added to the menu.
 */
void ConfigEditorController::fillInheritanceMenu(QMenu& menu, const std::vector<Model::ConfigNode*>& list)
{
	// Cleanup all actions owned by the menu, to avoid that they accumulate
	// until the menu is destroyed (usually only at the end of the application).
	foreach(QAction* action, menu.actions())
	{
		delete action;
	}
	menu.clear();

	// fills the menu with an action for each element in the list.
	for(std::vector<Model::ConfigNode*>::const_iterator iNode = list.begin();
		iNode != list.end(); ++iNode)
	{
		std::string nodePath = (*iNode)->getPath().toString();
		QAction* action = new QAction(QString::fromStdString(nodePath), &menu);
		menu.addAction(action);
		connect(action, SIGNAL(triggered()), this, SLOT(onNodeRefClicked()));
	}

	menu.setEnabled(list.size() > 0);
}

/**
 * Returns true if a selection has more than 1 item.
 */
bool ConfigEditorController::isMultiselection(const QItemSelection& selection) const
{
	int selectedCount = 0;
	for(QItemSelection::const_iterator iRange = selection.begin(); iRange != selection.end(); ++iRange)
	{
		for(int row = iRange->top(); row <= iRange->bottom(); ++row)
		{
			++selectedCount;
			if(selectedCount > 1)
			{
				return true;
			}
		}
	}
	return false;
}

}
}
