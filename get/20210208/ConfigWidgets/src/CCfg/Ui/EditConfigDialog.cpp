/*
 * $Id: EditConfigDialog.cpp 1585 2018-11-12 16:08:39Z psizun $
 * @file EditConfigDialog.cpp
 * @created 3 mai 2012
 * @author sizun
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * -----------------------------------------------------------------------------
 */

#include "CCfg/Ui/EditConfigDialog.h"
#include "CCfg/Ui/EditConfigDialog.ui.h"

#include "CCfg/Ui/ErrorGuiPrinter.h"
#include "CCfg/Ui/AttributeBuilderWidget.h"
#include "CCfg/Ui/ConfigEditorTree.h"
#include "CCfg/Ui/ClosableTabWidget.h"
#include "CCfg/Ui/AttributeFinderDialog.h"
#include "CCfg/Control/Workspace.h"
#include "CCfg/Control/BaseConfigEditorController.h"
#include "CCfg/Io/Url.h"
#include "CCfg/Model/ConfigModel.h"
#include "Utils/Error/AbstractPrinter.h"
#include "Utils/Error/ExcepHandling.h"

#include <QAction>
#include <QBoxLayout>
#include <QCloseEvent>
#include <QDir>
#include <QFileInfo>
#include <QIcon>
#include <QLayout>
#include <QMessageBox>
#include <QStatusBar>
#include <QToolButton>
#include <QUndoStack>
#include <QVBoxLayout>
#include <memory>

namespace CCfg {
namespace Ui {
using CCfg::Model::ConfigModel;
using namespace Utils::Error;
//_____________________________________________________________________________
/**
 * Default constructor.
 */
EditConfigDialog::EditConfigDialog(QWidget* parentWidget, Qt::WindowFlags flags)
	: QDialog(parentWidget, flags), widgets(new ::Ui::EditConfigDialog), workspaceCtrl(0)
{
	Reporter::setPrinter(std::auto_ptr<AbstractPrinter>(new ErrorGuiPrinter(this)));
	widgets->setupUi(this);
	statusBar = new QStatusBar(this);
	QBoxLayout* boxLayout = qobject_cast<QBoxLayout*>(layout());
	boxLayout->insertWidget(boxLayout->indexOf(widgets->buttonBox), statusBar);
	createActions();
	connectActions();
	setupUndoFramework();
	initAttributeBuilder();
	initAttributeFinder();
	initToolBar();

	// Prevent user from closing config
	QToolButton* closeButton = qobject_cast<QToolButton *>(getWorkspaceUi()->cornerWidget());
	closeButton->setDisabled(true);
	closeButton->setAutoRaise(false);
	closeButton->defaultAction()->disconnect();
}
//_____________________________________________________________________________
/**
 * Destructor.
 */
EditConfigDialog::~EditConfigDialog()
{
}
//_____________________________________________________________________________
void EditConfigDialog::reject()
{
	try
	{
		if(workspaceCtrl == 0 or getWorkspace()->closeAll())
		{
			QDialog::reject();
			return;
		}
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
/**
 * Returns whether all configurations in the workspace are clean.
 */
bool EditConfigDialog::isWorkspaceClean() const
{
	bool isClean = true;
	foreach(QUndoStack* stack, undoGroup.stacks())
	{
		isClean = isClean and stack->isClean();
	}
	return isClean;
}
//_____________________________________________________________________________

/**
 * Returns the workspace controller (creating it if necessary).
 */
Control::Workspace* EditConfigDialog::getWorkspace()
{
	if(workspaceCtrl == 0)
	{
		initWorkspace();
	}
	return workspaceCtrl;
}
//_____________________________________________________________________________

/**
 * Returns the workspace user interface.
 */
Ui::ClosableTabWidget* EditConfigDialog::getWorkspaceUi()
{
	return widgets->workspaceUi;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionOpen() const
{
	return actionOpen;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionSave() const
{
	return actionSave;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionSaveAs() const
{
	return actionSaveAs;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionSaveAll() const
{
	return actionSaveAll;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionShowName() const
{
	return actionShowName;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionShowValue() const
{
	return actionShowValue;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionShowUnit() const
{
	return actionShowUnit;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionShowRange() const
{
	return actionShowRange;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionShowInclude() const
{
	return actionShowInclude;
}
//_____________________________________________________________________________
QAction* EditConfigDialog::getActionShowBacktrackLevel() const
{
	return actionShowBacktrackLevel;
}
//_____________________________________________________________________________
/**
 * Create actions.
 */
void EditConfigDialog::createActions()
{
	actionNew = new QAction(QIcon(":/oxygen/document-new"), tr("New"), this);
	actionNew->setShortcuts(QKeySequence::New);

	actionOpen = new QAction(QIcon(":/oxygen/document-open"), tr("Open"), this);
	actionOpen->setShortcuts(QKeySequence::Open);

	actionSave = new QAction(QIcon(":/oxygen/document-save"), tr("Save"), this);
	actionSave->setShortcuts(QKeySequence::Save);

	actionSaveAs = new QAction(QIcon(":/oxygen/document-save-as"), tr("Save As"), this);
	actionSaveAs->setShortcuts(QKeySequence::SaveAs);

	actionSaveAll = new QAction(QIcon(":/oxygen/document-save-all"), tr("Save All"), this);

	actionShowValue = new QAction(tr("Show value"), this);
	actionShowUnit = new QAction(tr("Show unit"), this);
	actionShowRange = new QAction(tr("Show range"), this);
	actionShowBacktrackLevel = new QAction(tr("Show backtrack level"), this);
	actionShowInclude = new QAction(tr("Show include"), this);

	actionShowAttributeBuilder = new QAction(QIcon(":/builder/default_star.png"), tr("Show attibute builder"), this);
#if QT_VERSION >= 0x040600
	actionShowAttributeBuilder->setIcon(QIcon::fromTheme("applications-utilities", QIcon(":/builder/default_star.png")));
#endif
}
//_____________________________________________________________________________
/**
 * Connects menu actions signals to their slots.
 */
void EditConfigDialog::connectActions()
{
	connect(actionNew, SIGNAL(triggered()), this, SLOT(newConfig()));
	connect(actionOpen, SIGNAL(triggered()), this, SLOT(loadFile()));
	connect(actionSave, SIGNAL(triggered()), this, SLOT(save()));
	connect(actionSaveAs, SIGNAL(triggered()), this, SLOT(saveAs()));
	connect(actionSaveAll, SIGNAL(triggered()), this, SLOT(saveAll()));

	connect(actionShowValue, SIGNAL(triggered(bool)),
			this, SLOT(showValueColumn(bool)));
	connect(actionShowUnit, SIGNAL(triggered(bool)),
			this, SLOT(showUnitColumn(bool)));
	connect(actionShowRange, SIGNAL(triggered(bool)),
			this, SLOT(showRangeColumn(bool)));
	connect(actionShowBacktrackLevel, SIGNAL(triggered(bool)),
			this, SLOT(showBacktrackLevelColumn(bool)));
	connect(actionShowInclude, SIGNAL(triggered(bool)),
			this, SLOT(showIncludeUrlColumn(bool)));

	connect(actionShowAttributeBuilder, SIGNAL(triggered(bool)),
				this, SLOT(toggleAttributeBuilder()));
}
//_____________________________________________________________________________
/**
 * Creates tool bar.
 */
void EditConfigDialog::initToolBar()
{
	widgets->buttonUndo->setDefaultAction(actionUndo);
	widgets->buttonRedo->setDefaultAction(actionRedo);
	widgets->buttonBuilder->setDefaultAction(actionShowAttributeBuilder);
	widgets->buttonFind->setDefaultAction(actionFind);
}
//_____________________________________________________________________________
/**
 * Initializes the Attribute Builder container widget.
 */
void EditConfigDialog::initAttributeBuilder()
{
	widgets->attributeBuilder->parentWidget()->hide();
}
//_____________________________________________________________________________
/**
 * Intializes the Workspace (tabbed multidocument interface for editing configurations).
 * It must not be called in the constructor, otherwise the virtual method
 * createWorkspace won't be overriden.
 */
void EditConfigDialog::initWorkspace()
{
	//widgets->workspaceUi->setCornerWidget(0);

	workspaceCtrl = createWorkspace(widgets->workspaceUi);
	workspaceCtrl->setAttributeBuilder(widgets->attributeBuilder);

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
//_____________________________________________________________________________
/**
 * Initializes the Attribute Finder dialog.
 */
void EditConfigDialog::initAttributeFinder()
{
	attributeFinder = new AttributeFinderDialog(this);
	attributeFinder->setModal(false);
	attributeFinder->hide();

	actionFind = new QAction(tr("Find..."), this);
#if QT_VERSION >= 0x040600
	actionFind->setIcon(QIcon::fromTheme("edit-find"));
#else
	actionFind->setIcon(QIcon(":/oxygen/edit-find"))
#endif
	actionFind->setEnabled(false);
	addAction(actionFind);

	actionFind->setShortcuts(QKeySequence::Find);
	actionFindNext = new QAction(tr("Find Next"), this);
	actionFindNext->setShortcuts(QKeySequence::FindNext);
	actionFindNext->setEnabled(false);
	addAction(actionFindNext);

	connect(attributeFinder, SIGNAL(searchAttribute(const QString &, bool, bool, int)), this, SLOT(onSearchAttribute(const QString &, bool, bool, int)));
	connect(actionFind, SIGNAL(triggered()), this, SLOT(findAttribute()));
	connect(actionFindNext, SIGNAL(triggered()), this, SLOT(findNextAttribute()));
}
//_____________________________________________________________________________
/**
 * Initializes the undo framework to integrate it in the UI.
 */
void EditConfigDialog::setupUndoFramework()
{
	actionUndo = undoGroup.createUndoAction(this);
	actionRedo = undoGroup.createRedoAction(this);
	actionUndo->setIcon(QIcon(":/oxygen/edit-undo"));
	actionRedo->setIcon(QIcon(":/oxygen/edit-redo"));
	actionUndo->setShortcut(QKeySequence::Undo);
	actionRedo->setShortcut(QKeySequence::Redo);

	widgets->buttonBox->addAction(actionUndo);
	widgets->buttonBox->addAction(actionRedo);

    connect(&undoGroup, SIGNAL(cleanChanged(bool)), this, SLOT(onCleanChanged(bool)));
}
//_____________________________________________________________________________
/**
 * Creates a new configuration in the Workspace.
 */
void EditConfigDialog::newConfig()
{
	try
	{
		getWorkspace()->newConfig();
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
/**
 * Opens a dialog to edit the configuration file whose path is provided.
 */
void EditConfigDialog::editConfig(QWidget* parentWidget, const std::string & cfgFilePath, const QString & title)
{
	EditConfigDialog dialog(parentWidget);
	if (not title.isEmpty())
		dialog.setWindowTitle(title);
	dialog.openConfig(cfgFilePath);
	dialog.exec();
}
//_____________________________________________________________________________
/**
 * Loads a configuration file in the workspace.
 */
void EditConfigDialog::openConfig(const std::string & cfgFilePath)
{
	if (not cfgFilePath.empty())
	{
		getWorkspace()->openConfig(CCfg::Io::Url::fromPath(cfgFilePath))->getModel();
	}
}
//_____________________________________________________________________________
/**
 * Loads the configuration files selected by the user in the Workspace.
 */
void EditConfigDialog::loadFile()
{
	try
	{
		getWorkspace()->open();
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
/**
 * Saves the configuration in the active tab.
 */
bool EditConfigDialog::save()
{
	try
	{
		return getWorkspace()->save();
	}
	CATCH_AND_REPORT
	return false;
}
//_____________________________________________________________________________
/**
 * Saves the active configuration in another file.
 */
bool EditConfigDialog::saveAs()
{
	try
	{
		return getWorkspace()->saveAs();
	}
	CATCH_AND_REPORT
	return false;
}
//_____________________________________________________________________________
/**
 * Saves all configuration.
 */
void EditConfigDialog::saveAll()
{
	try
	{
		return getWorkspace()->saveAll();
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
/**
 * Shows or hide the attribute builder.
 */
void EditConfigDialog::toggleAttributeBuilder()
{
	bool isVisible = widgets->attributeBuilderParentWidget->isVisible();
	widgets->attributeBuilderParentWidget->setVisible(not isVisible);
}
//_____________________________________________________________________________
/**
 * Shows or hide the specified column of the current tab configuration tree.
 */
void EditConfigDialog::showColumn(int index, bool show)
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
//_____________________________________________________________________________
/**
 * Shows or hide the value column of the current tab configuration tree.
 */
void EditConfigDialog::showValueColumn(bool show)
{
	showColumn(ConfigModel::ValueColumn, show);
}
//_____________________________________________________________________________
/**
 * Shows or hide the unit column of the current tab configuration tree.
 */
void EditConfigDialog::showUnitColumn(bool show)
{
	showColumn(ConfigModel::UnitColumn, show);
}
//_____________________________________________________________________________
/**
 * Shows or hide the range column of the current tab configuration tree.
 */
void EditConfigDialog::showRangeColumn(bool show)
{
	showColumn(ConfigModel::RangeColumn, show);
}
//_____________________________________________________________________________
/**
 * Shows or hide the backtrack level column of the current tab configuration tree.
 */
void EditConfigDialog::showBacktrackLevelColumn(bool show)
{
	showColumn(ConfigModel::BacktrackColumn, show);
}
//_____________________________________________________________________________
/**
 * Shows or hide the include url column of the current tab configuration tree.
 */
void EditConfigDialog::showIncludeUrlColumn(bool show)
{
	showColumn(ConfigModel::IncludeColumn, show);
}
//_____________________________________________________________________________
/**
 * Displays the attribute finder dialog.
 */
void EditConfigDialog::findAttribute()
{
	attributeFinder->show();
}
//_____________________________________________________________________________
/**
 * Repeats last search.
 */
void EditConfigDialog::findNextAttribute()
{
	attributeFinder->findAttribute();
}
//_____________________________________________________________________________
/**
 * Searches for an attribute matching given criteria.
 */
void EditConfigDialog::onSearchAttribute(const QString & searchText, bool caseSens, bool exactMatch, int column)
{
	try
	{
		ConfigEditorTree* tree = getWorkspace()->getCurrentController()->getView();
		if (not tree) return;
		tree->findAttribute(searchText, caseSens, exactMatch, column);
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
/**
 * Creates a new instance of workspace controller.
 * @param ui the user interface of the workspace.
 */
Control::Workspace* EditConfigDialog::createWorkspace(ClosableTabWidget* ui)
{
	return new Control::Workspace(ui, this);
}
//_____________________________________________________________________________
/**
 * Reacts to the event closing this main window.
 * It asks the user to save its documents. The user can interrupt the operation
 * if its unsaved documents requires additional processing before saving.
 */
void EditConfigDialog::closeEvent(QCloseEvent* event)
{
	try
	{
		if(workspaceCtrl == 0 or getWorkspace()->closeAll())
		{
			event->accept();
			return;
		}
	}
	CATCH_AND_REPORT
	event->ignore();
}
//_____________________________________________________________________________
/**
 * Reacts to the insertion of a new model in the Workspace.
 */
void EditConfigDialog::onModelAdded(Model::ConfigModel* model)
{
	try
	{
		undoGroup.addStack(model->getUndoStack());
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
/**
 * Reacts to the removal of a model in the Workspace.
 * A model is removed when no tab display it any more.
 */
void EditConfigDialog::onModelRemoved(Model::ConfigModel* model)
{
	try
	{
		undoGroup.removeStack(model->getUndoStack());
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
/**
 * Reacts to the activation of a model in the Workspace.
 * Note that the signal is not emitted when a tab diplaying the same model as
 * the previous is activated.
 */
void EditConfigDialog::onModelActivated(Model::ConfigModel* current, Model::ConfigModel* /*previous*/)
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
//_____________________________________________________________________________
/**
 * Reacts to a change of activated tab.
 */
void EditConfigDialog::onTabActivated(Control::BaseConfigEditorController* current, Control::BaseConfigEditorController* /*previous*/)
{
	try
	{
		bool hasCurrent = (current != 0);
		if(hasCurrent)
		{
			ConfigEditorTree* tree = current->getView();
			actionShowValue->setChecked(not tree->isColumnHidden(ConfigModel::ValueColumn));
			actionShowUnit->setChecked(not tree->isColumnHidden(ConfigModel::UnitColumn));
			actionShowRange->setChecked(not tree->isColumnHidden(ConfigModel::RangeColumn));
			actionShowBacktrackLevel->setChecked(not tree->isColumnHidden(ConfigModel::BacktrackColumn));
			actionShowInclude->setChecked(not tree->isColumnHidden(ConfigModel::IncludeColumn));
		}
		actionSaveAs->setEnabled(hasCurrent);
		actionFind->setEnabled(hasCurrent);
		actionFindNext->setEnabled(hasCurrent);
		actionShowValue->setEnabled(hasCurrent);
		actionShowUnit->setEnabled(hasCurrent);
		actionShowRange->setEnabled(hasCurrent);
		actionShowBacktrackLevel->setEnabled(hasCurrent);
		actionShowInclude->setEnabled(hasCurrent);
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
/**
 * Reacts to a change of the clean state of the undo group
 * It means that the active model was modified or that a dirty model has been
 * activated.
 */
void EditConfigDialog::onCleanChanged(bool isClean)
{
	try
	{
		actionSave->setEnabled(not isClean);
		actionSaveAll->setEnabled(not isWorkspaceClean());
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
/**
 * Updates the status bar when the current node changes.
 */
void EditConfigDialog::onCurrentNodeChanged(Model::ConfigNode* current,
		Model::ConfigModel* /*currentModel*/,
		Model::ConfigNode* /*previous*/, Model::ConfigModel* /*previousModel*/)
{
	try
	{
		if(current != 0)
		{
			statusBar->showMessage(QString::fromStdString(current->getPath().toString()), 0);
		}
	}
	CATCH_AND_REPORT
}
//_____________________________________________________________________________
} /* namespace Ui */
} /* namespace CCfg */
