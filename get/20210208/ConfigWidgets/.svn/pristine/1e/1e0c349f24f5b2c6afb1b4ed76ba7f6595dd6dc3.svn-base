#include "EditorMainWindow.h"
#include "ErrorGuiPrinter.h"
#include "../Control/Workspace.h"
#include "../Control/BaseConfigEditorController.h"
#include "AttributeBuilderWidget.h"
#include "AttributeFinderDialog.h"
#include "ConfigEditorTree.h"
#include "ClosableTabWidget.h"
#include "../Model/ConfigModel.h"
#include "Utils/Error/ExcepHandling.h"
#include "CCfg/Ui/MainWindow.ui.h"
#include "CCfg/Io/Url.h"
#include <utl/Logging.h>
#include <QDockWidget>
#include <QLayout>
#include <QCloseEvent>
#include <QUndoStack>
#include <QMessageBox>
#include <QSettings>
#include <QAction>
#include <QDir>
#include <QIcon>

namespace CCfg
{
namespace Ui
{
using CCfg::Model::ConfigModel;
using namespace Utils::Error;

///////////////////////////////////////////////////////////////////////////////
// EditorMainWindow class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
EditorMainWindow::EditorMainWindow(QWidget* parentWidget)
	: Utils::Ui::DesignedWidget<Ui_MainWindow, QMainWindow>(parentWidget),
	workspaceCtrl(0), workspaceUi(0), attributeBuilder(0),
	attributeFinder(0), undoGroup(), actionUndo(0), actionRedo(0)
{
	Reporter::setPrinter(std::auto_ptr<AbstractPrinter>(new ErrorGuiPrinter(this)));
//#if QT_VERSION >= 0x040600
//	ui_->actionFind->setIcon(QIcon::fromTheme("edit-find"));
//#endif
	ui_->actionAbout->setShortcuts(QKeySequence::HelpContents);
#if QT_VERSION >= 0x040600
	ui_->actionAbout->setIcon(QIcon::fromTheme("help-about"));
#endif
	ui_->actionNew->setShortcuts(QKeySequence::New);
	ui_->actionOpen->setShortcuts(QKeySequence::Open);
	ui_->actionSave->setShortcuts(QKeySequence::Save);
	ui_->actionSaveAs->setShortcuts(QKeySequence::SaveAs);
	ui_->actionQuit->setShortcuts(QKeySequence::Quit);
	ui_->actionFind->setShortcuts(QKeySequence::Find);
	ui_->actionFindNext->setShortcuts(QKeySequence::FindNext);

	restoreWindowShape();
	connectActions();
	setupUndoFramework();
}

/**
 * Destructor.
 */
EditorMainWindow::~EditorMainWindow()
{
}

void EditorMainWindow::setAttributeBuilder(AttributeBuilderWidget* attributeBuilder)
{
	this->attributeBuilder = attributeBuilder;
}

/**
 * Returns whether all configurations in the workspace are clean.
 */
bool EditorMainWindow::isWorkspaceClean() const
{
	bool isClean = true;
	foreach(QUndoStack* stack, undoGroup.stacks())
	{
		isClean = isClean and stack->isClean();
	}
	return isClean;
}

/**
 * Returns the workspace controller (creating it if necessary).
 */
Control::Workspace* EditorMainWindow::getWorkspace()
{
	if(workspaceCtrl == 0)
	{
		initWorkspace();
	}
	return workspaceCtrl;
}

AttributeBuilderWidget* EditorMainWindow::getAttributeBuilder()
{
	if(attributeBuilder == 0)
	{
		initAttributeBuilder();
	}
	return attributeBuilder;
}

AttributeFinderDialog* EditorMainWindow::getAttributeFinder()
{
	if(attributeFinder == 0)
	{
		initAttributeFinder();
	}
	return attributeFinder;
}

/**
 * Returns the workspace user interface (creating it if necessary).
 */
Ui::ClosableTabWidget* EditorMainWindow::getWorkspaceUi()
{
	if(workspaceUi == 0)
	{
		initWorkspace();
	}
	return workspaceUi;
}

QMenuBar* EditorMainWindow::getMenubar() const
{
	return ui_->menubar;
}

QStatusBar* EditorMainWindow::getStatusbar() const
{
	return ui_->statusbar;
}

QToolBar* EditorMainWindow::getToolBar() const
{
	return ui_->toolBar;
}

QMenu* EditorMainWindow::getMenuFile() const
{
	return ui_->menuFile;
}

QMenu* EditorMainWindow::getMenuEdit() const
{
	return ui_->menuEdit;
}

QMenu* EditorMainWindow::getMenuView() const
{
	return ui_->menuView;
}

QMenu* EditorMainWindow::getMenuHelp() const
{
	return ui_->menuHelp;
}

QMenu* EditorMainWindow::getMenuColumns() const
{
	return ui_->menuColumns;
}

QAction* EditorMainWindow::getActionNew() const
{
	return ui_->actionNew;
}

QAction* EditorMainWindow::getActionOpen() const
{
	return ui_->actionOpen;
}

QAction* EditorMainWindow::getActionSave() const
{
	return ui_->actionSave;
}

QAction* EditorMainWindow::getActionSaveAs() const
{
	return ui_->actionSaveAs;
}

QAction* EditorMainWindow::getActionSaveAll() const
{
	return ui_->actionSaveAll;
}

QAction* EditorMainWindow::getActionQuit() const
{
	return ui_->actionQuit;
}

QAction* EditorMainWindow::getActionFind() const
{
	return ui_->actionFind;
}

QAction* EditorMainWindow::getActionShowName() const
{
	return ui_->actionShowName;
}

QAction* EditorMainWindow::getActionShowValue() const
{
	return ui_->actionShowValue;
}

QAction* EditorMainWindow::getActionShowUnit() const
{
	return ui_->actionShowUnit;
}

QAction* EditorMainWindow::getActionShowRange() const
{
	return ui_->actionShowRange;
}

QAction* EditorMainWindow::getActionShowInclude() const
{
	return ui_->actionShowInclude;
}

QAction* EditorMainWindow::getActionShowBacktrackLevel() const
{
	return ui_->actionShowBacktrackLevel;
}

QAction* EditorMainWindow::getActionAbout() const
{
	return ui_->actionAbout;
}

/**
 * Connects menu actions signals to their slots.
 */
void EditorMainWindow::connectActions()
{
	connect(ui_->actionNew, SIGNAL(triggered()), this, SLOT(newConfig()));
	connect(ui_->actionOpen, SIGNAL(triggered()), this, SLOT(openFile()));
	connect(ui_->actionSave, SIGNAL(triggered()), this, SLOT(save()));
	connect(ui_->actionSaveAs, SIGNAL(triggered()), this, SLOT(saveAs()));
	connect(ui_->actionSaveAll, SIGNAL(triggered()), this, SLOT(saveAll()));
	connect(ui_->actionQuit, SIGNAL(triggered()), qApp, SLOT(closeAllWindows()));
	connect(ui_->actionFind, SIGNAL(triggered()), this, SLOT(findAttribute()));
	connect(ui_->actionFindNext, SIGNAL(triggered()), this, SLOT(findNextAttribute()));
	connect(ui_->actionAbout, SIGNAL(triggered()), this, SLOT(about()));

	connect(ui_->actionShowValue, SIGNAL(triggered(bool)),
			this, SLOT(showValueColumn(bool)));
	connect(ui_->actionShowUnit, SIGNAL(triggered(bool)),
			this, SLOT(showUnitColumn(bool)));
	connect(ui_->actionShowRange, SIGNAL(triggered(bool)),
			this, SLOT(showRangeColumn(bool)));
	connect(ui_->actionShowBacktrackLevel, SIGNAL(triggered(bool)),
			this, SLOT(showBacktrackLevelColumn(bool)));
	connect(ui_->actionShowInclude, SIGNAL(triggered(bool)),
			this, SLOT(showIncludeUrlColumn(bool)));
}

/**
 * Initializes the Attribute Builder dock widget.
 */
void EditorMainWindow::initAttributeBuilder()
{
	QDockWidget* attributeBuilderDock = new QDockWidget(tr("Attribute builder"), this);
	attributeBuilder = new AttributeBuilderWidget(attributeBuilderDock);
	attributeBuilderDock->setObjectName("attributeBuilderDock");
	attributeBuilderDock->setWidget(attributeBuilder);
	attributeBuilderDock->layout()->setAlignment(attributeBuilder, Qt::AlignTop);
	attributeBuilderDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	addDockWidget(Qt::LeftDockWidgetArea, attributeBuilderDock, Qt::Vertical);
	attributeBuilderDock->hide();

	ui_->menuView->addAction(attributeBuilderDock->toggleViewAction());
}

/**
 * Initializes the Attribute Finder dialog.
 */
void EditorMainWindow::initAttributeFinder()
{
	attributeFinder = new AttributeFinderDialog(this);
	attributeFinder->setModal(false);
	attributeFinder->hide();

	connect(attributeFinder, SIGNAL(searchAttribute(const QString &, bool, bool, int)), this, SLOT(onSearchAttribute(const QString &, bool, bool, int)));
}

/**
 * Initializes the Workspace (tabbed multidocument interface for editing configurations).
 * It must not be called in the constructor, otherwise the virtual method
 * createWorkspace won't be overridden.
 */
void EditorMainWindow::initWorkspace()
{
	workspaceUi = new ClosableTabWidget(this);
	setCentralWidget(workspaceUi);

	workspaceCtrl = createWorkspace(workspaceUi);
	workspaceCtrl->setAttributeBuilder(getAttributeBuilder());

	connect(workspaceCtrl, SIGNAL(modelAdded(Model::ConfigModel*)),
			this, SLOT(onModelAdded(Model::ConfigModel*)));
	connect(workspaceCtrl, SIGNAL(modelRemoved(Model::ConfigModel*)),
			this, SLOT(onModelRemoved(Model::ConfigModel*)));
	connect(workspaceCtrl, SIGNAL(modelActivated(Model::ConfigModel*, Model::ConfigModel*)),
			this, SLOT(onModelActivated(Model::ConfigModel*, Model::ConfigModel*)));
	connect(workspaceCtrl, SIGNAL(currentNodeChanged(Model::ConfigNode*, Model::ConfigModel*, Model::ConfigNode*, Model::ConfigModel*)),
			this, SLOT(onCurrentNodeChanged(Model::ConfigNode*, Model::ConfigModel*, Model::ConfigNode*, Model::ConfigModel*)));
	connect(workspaceCtrl, SIGNAL(tabActivated(Control::BaseConfigEditorController*, Control::BaseConfigEditorController*)),
			this, SLOT(onTabActivated(Control::BaseConfigEditorController*, Control::BaseConfigEditorController*)));
}

/**
 * Initializes the undo framework to integrate it in the UI.
 */
void EditorMainWindow::setupUndoFramework()
{
	actionUndo = undoGroup.createUndoAction(this);
	actionRedo = undoGroup.createRedoAction(this);
	actionUndo->setIcon(QIcon(":/oxygen/edit-undo"));
	actionRedo->setIcon(QIcon(":/oxygen/edit-redo"));
	actionUndo->setShortcut(QKeySequence::Undo);
	actionRedo->setShortcut(QKeySequence::Redo);

	ui_->menuEdit->addSeparator();
	ui_->menuEdit->addAction(actionUndo);
	ui_->menuEdit->addAction(actionRedo);

	ui_->toolBar->addSeparator();
    ui_->toolBar->addAction(actionUndo);
    ui_->toolBar->addAction(actionRedo);

    connect(&undoGroup, SIGNAL(cleanChanged(bool)), this, SLOT(onCleanChanged(bool)));
}

/**
 * Creates a new configuration in the Workspace.
 */
void EditorMainWindow::newConfig()
{
	try
	{
		getWorkspace()->newConfig();
	}
	CATCH_AND_REPORT
}

/**
 * Opens the configuration files selected by the user in the Workspace.
 */
void EditorMainWindow::openFile()
{
	try
	{
		getWorkspace()->open();
	}
	CATCH_AND_REPORT
}

/**
 * Opens the specified configuration files.
 */
void EditorMainWindow::openFiles(const QStringList& fileList)
{
	try
	{
		getWorkspace()->openFiles(fileList);
	}
	CATCH_AND_REPORT
}

/**
 * Saves the configuration in the active tab.
 */
bool EditorMainWindow::save()
{
	try
	{
		return getWorkspace()->save();
	}
	CATCH_AND_REPORT
	return false;
}

/**
 * Saves the active configuration in another file.
 */
bool EditorMainWindow::saveAs()
{
	try
	{
		return getWorkspace()->saveAs();
	}
	CATCH_AND_REPORT
	return false;
}

/**
 * Saves all configuration.
 */
void EditorMainWindow::saveAll()
{
	try
	{
		return getWorkspace()->saveAll();
	}
	CATCH_AND_REPORT
}

/**
 * Shows or hide the specified column of the current tab configuration tree.
 */
void EditorMainWindow::showColumn(int index, bool show)
{
	try
	{
		ConfigEditorTree* tree = getWorkspace()->getCurrentController()->getView();
		if(tree != 0)
		{
			tree->setColumnHidden(index, not show);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Shows or hide the value column of the current tab configuration tree.
 */
void EditorMainWindow::showValueColumn(bool show)
{
	showColumn(ConfigModel::ValueColumn, show);
}

/**
 * Shows or hide the unit column of the current tab configuration tree.
 */
void EditorMainWindow::showUnitColumn(bool show)
{
	showColumn(ConfigModel::UnitColumn, show);
}

/**
 * Shows or hide the range column of the current tab configuration tree.
 */
void EditorMainWindow::showRangeColumn(bool show)
{
	showColumn(ConfigModel::RangeColumn, show);
}

/**
 * Shows or hide the backtrack level column of the current tab configuration tree.
 */
void EditorMainWindow::showBacktrackLevelColumn(bool show)
{
	showColumn(ConfigModel::BacktrackColumn, show);
}

/**
 * Shows or hide the include url column of the current tab configuration tree.
 */
void EditorMainWindow::showIncludeUrlColumn(bool show)
{
	showColumn(ConfigModel::IncludeColumn, show);
}

/**
 * Creates a new instance of workspace controller.
 * @param ui the user interface of the workspace.
 */
Control::Workspace* EditorMainWindow::createWorkspace(ClosableTabWidget* ui)
{
	return new Control::Workspace(ui, this);
}

void EditorMainWindow::saveWindowShape()
{
	LOG_DEBUG() << "saving EditorMainWindow shape";
	QSettings settings;
	settings.beginGroup("windowShape");
	settings.beginGroup("EditorMainWindow");
    settings.setValue("geometry", saveGeometry());
    settings.setValue("state", saveState());
    settings.endGroup();
    settings.endGroup();
}

void EditorMainWindow::restoreWindowShape()
{
	LOG_DEBUG() << "restoring EditorMainWindow shape";
	QSettings settings;
	settings.beginGroup("windowShape");
	settings.beginGroup("EditorMainWindow");
	restoreGeometry(settings.value("geometry").toByteArray());
	restoreState(settings.value("state").toByteArray());
    settings.endGroup();
    settings.endGroup();
}

/**
 * Reacts to the event closing this main window.
 * It asks the user to save its documents. The user can interrupt the operation
 * if its unsaved documents requires additional processing before saving.
 */
void EditorMainWindow::closeEvent(QCloseEvent* event)
{
	try
	{
		if(workspaceCtrl == 0 or getWorkspace()->closeAll())
		{
			saveWindowShape();
			event->accept();
			return;
		}
	}
	CATCH_AND_REPORT
	event->ignore();
}

/**
 * Reacts to the insertion of a new model in the Workspace.
 */
void EditorMainWindow::onModelAdded(Model::ConfigModel* model)
{
	try
	{
		undoGroup.addStack(model->getUndoStack());
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the removal of a model in the Workspace.
 * A model is removed when no tab display it any more.
 */
void EditorMainWindow::onModelRemoved(Model::ConfigModel* model)
{
	try
	{
		undoGroup.removeStack(model->getUndoStack());
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the activation of a model in the Workspace.
 * Note that the signal is not emitted when a tab diplaying the same model as
 * the previous is activated.
 */
void EditorMainWindow::onModelActivated(Model::ConfigModel* current, Model::ConfigModel* /*previous*/)
{
	try
	{
		if(current != 0)
		{
			undoGroup.setActiveStack(current->getUndoStack());
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a change of activated tab.
 */
void EditorMainWindow::onTabActivated(Control::BaseConfigEditorController* current, Control::BaseConfigEditorController* /*previous*/)
{
	try
	{
		bool hasCurrent = (current != 0);
		if(hasCurrent)
		{
			ConfigEditorTree* tree = current->getView();
			ui_->actionShowValue->setChecked(not tree->isColumnHidden(ConfigModel::ValueColumn));
			ui_->actionShowUnit->setChecked(not tree->isColumnHidden(ConfigModel::UnitColumn));
			ui_->actionShowRange->setChecked(not tree->isColumnHidden(ConfigModel::RangeColumn));
			ui_->actionShowBacktrackLevel->setChecked(not tree->isColumnHidden(ConfigModel::BacktrackColumn));
			ui_->actionShowInclude->setChecked(not tree->isColumnHidden(ConfigModel::IncludeColumn));
		}
		else
		{
			if(attributeFinder != 0)
			{
				attributeFinder->hide();
			}
		}
		ui_->actionSaveAs->setEnabled(hasCurrent);
		ui_->actionFind->setEnabled(hasCurrent);
		ui_->actionFindNext->setEnabled(hasCurrent);
		ui_->actionShowValue->setEnabled(hasCurrent);
		ui_->actionShowUnit->setEnabled(hasCurrent);
		ui_->actionShowRange->setEnabled(hasCurrent);
		ui_->actionShowBacktrackLevel->setEnabled(hasCurrent);
		ui_->actionShowInclude->setEnabled(hasCurrent);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a change of the clean state of the undo group
 * It means that the active model was modified or that a dirty model has been
 * activated.
 */
void EditorMainWindow::onCleanChanged(bool isClean)
{
	try
	{
		ui_->actionSave->setEnabled(not isClean);
		ui_->actionSaveAll->setEnabled(not isWorkspaceClean());
	}
	CATCH_AND_REPORT
}

/**
 * Updates the status bar when the current node changes.
 */
void EditorMainWindow::onCurrentNodeChanged(Model::ConfigNode* current,
		Model::ConfigModel* /*currentModel*/,
		Model::ConfigNode* /*previous*/, Model::ConfigModel* /*previousModel*/)
{
	try
	{
		if(current != 0)
		{
			ui_->statusbar->showMessage(QString::fromStdString(current->getPath().toString()), 0);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Displays the attribute finder dialog.
 */
void EditorMainWindow::findAttribute()
{
	getAttributeFinder()->show();
}

/**
 * Repeats last search.
 */
void EditorMainWindow::findNextAttribute()
{
	getAttributeFinder()->findAttribute();
}

/**
 * Searches for an attribute matching given criteria.
 * @param searchText String to search for.
 * @param caseSens Whether the search should be case sensitive or not.
 * @param exactMatch Whether the match should simply contain the pattern or be an exact match.
 * @param column Index of the column to search in.
 */
void EditorMainWindow::onSearchAttribute(const QString & searchText, bool caseSens, bool exactMatch, int column)
{
	try
	{
		Control::BaseConfigEditorController* controller = getWorkspace()->getCurrentController();
		if (not controller) return;
		ConfigEditorTree* tree = controller->getView();
		if (not tree) return;
		tree->findAttribute(searchText, caseSens, exactMatch, column);
	}
	CATCH_AND_REPORT
}

/**
 * Displays a short information message about the configuration editor.
 */
void EditorMainWindow::about()
{
	QMessageBox::about(this, tr("About Configuration Editor"),
			tr("<h2>Compound Configuration Editor</h2>Copyright &copy; 2006-2012, <i>CEA, IRFU</i>"));
}

}
}
