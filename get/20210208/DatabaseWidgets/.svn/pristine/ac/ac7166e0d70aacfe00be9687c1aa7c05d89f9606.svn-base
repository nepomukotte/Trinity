#include "EditorMainWindow2.h"
#include "../Control/Workspace2.h"
#include "Utils/Error/ExcepHandling.h"
#include <QMenu>
#include <QAction>
#include <QMenuBar>

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::EditorMainWindow2 class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
EditorMainWindow2::EditorMainWindow2(QWidget* parentWidget) :
	EditorMainWindow(parentWidget),
	actionCheckout(0),
	actionCommit(0),
	actionCommitInto(0),
	actionRunCheckout(0),
	actionRunEdit(0)
{
	actionCheckout = new QAction(tr("Checkout..."), this);
	actionCheckout->setObjectName("actionCheckout");
	actionCheckout->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_C));

	actionCommit = new QAction(tr("Commit"), this);
	actionCommit->setObjectName("actionCommit");
	actionCommit->setEnabled(false);

	actionCommitInto = new QAction(tr("Commit into..."), this);
	actionCommitInto->setObjectName("actionCommitInto");
	actionCommitInto->setEnabled(false);

	actionRunCheckout = new QAction(tr("Checkout run..."), this);
	actionRunCheckout->setObjectName("actionRunCheckout");

	actionRunEdit = new QAction(tr("Edit run..."), this);
	actionRunEdit->setObjectName("actionRunEdit");

	QMenu* menuFile = getMenuFile();
	QAction* actionQuit = getActionQuit();
	menuFile->insertAction(actionQuit, actionCheckout);
	menuFile->insertAction(actionQuit, actionCommit);
	menuFile->insertAction(actionQuit, actionCommitInto);
	menuFile->insertSeparator(actionQuit);

	connect(actionCheckout, SIGNAL(triggered()), this, SLOT(checkout()));
	connect(actionCommit, SIGNAL(triggered()), this, SLOT(commit()));
	connect(actionCommitInto, SIGNAL(triggered()), this, SLOT(commitInto()));
	connect(actionRunCheckout, SIGNAL(triggered()), this, SLOT(checkoutRun()));
	connect(actionRunEdit, SIGNAL(triggered()), this, SLOT(editRun()));

	createMenuRun();
}

/**
 * Destructor.
 */
EditorMainWindow2::~EditorMainWindow2()
{
}

/**
 * Creates an additional menu for run management.
 */
void EditorMainWindow2::createMenuRun()
{
	// Insert menu "Run" before menu "Help"
	menuRun = new QMenu(tr("&Run"), this);
	menuBar()->insertMenu(getMenuHelp()->menuAction(), menuRun);

	menuRun->addAction(actionRunCheckout);
	menuRun->addAction(actionRunEdit);
}

/**
 * Returns the workspace controller (creating it if necessary).
 */
Control::Workspace2* EditorMainWindow2::getWorkspace()
{
	return dynamic_cast<Control::Workspace2*>(EditorMainWindow::getWorkspace());
}


/**
 * Displays a dialog box to select a configuration in the database, and checkout it.
 */
void EditorMainWindow2::checkout()
{
	try
	{
		getWorkspace()->checkoutConfig();
	}
	CATCH_AND_REPORT
}

/**
 * Commits the current configuration into the database from where it has been checked out.
 */
void EditorMainWindow2::commit()
{
	try
	{
		getWorkspace()->commit();
	}
	CATCH_AND_REPORT
}

/**
 * Displays a dialog box to select a configuration in the database where to commit
 * the current configuration.
 */
void EditorMainWindow2::commitInto()
{
	try
	{
		getWorkspace()->commitInto();
	}
	CATCH_AND_REPORT
}

/**
 * Displays a dialog box to checkout one or several run configurations.
 */
void EditorMainWindow2::checkoutRun()
{
	try
	{
		getWorkspace()->checkoutRun();
	}
	CATCH_AND_REPORT
}

/**
 * Displays a dialog box to add a configuration in the database to a run in the database.
 */
void EditorMainWindow2::editRun()
{
	try
	{
		getWorkspace()->editRun();
	}
	CATCH_AND_REPORT
}

/*
 * (see inherited doc)
 */
Control::Workspace* EditorMainWindow2::createWorkspace(ClosableTabWidget* ui)
{
	return new Control::Workspace2(ui, this);
}

/*
 * (see inherited doc)
 */
void EditorMainWindow2::onTabActivated(Control::BaseConfigEditorController* current, Control::BaseConfigEditorController* previous)
{
	try
	{
		bool hasCurrent = (current != 0);
		actionCommit->setEnabled(hasCurrent);
		actionCommitInto->setEnabled(hasCurrent);
	}
	CATCH_AND_REPORT
	EditorMainWindow::onTabActivated(current, previous);
}

}
}
