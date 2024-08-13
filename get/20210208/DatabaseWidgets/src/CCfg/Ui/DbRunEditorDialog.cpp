#include "DbRunEditorDialog.h"
#include "CCfg/Ui/DbRunEditorDialog.ui.h"
#include "CCfg/Model/AbstractDatabaseItem.h"
#include "CCfg/Model/RunItem.h"
#include "CCfg/Model/VersionItem.h"
#include "CCfg/Model/DatabaseModel.h"

namespace CCfg
{
using namespace Model;
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::DbRunEditorDialog class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
DbRunEditorDialog::DbRunEditorDialog(QWidget *parent)
    : QDialog(parent), widgets(new Ui_DbRunEditorDialog)
{
	widgets->setupUi(this);
	widgets->configBrowser->setFilterMode(DatabaseFilterProxyModel::ConfigOnly);
	widgets->configBrowser->setAcceptMode(DatabaseBrowserTree::AcceptLoad);
	widgets->runBrowser->setFilterMode(DatabaseFilterProxyModel::RunOnly);
	widgets->runBrowser->setAcceptMode(DatabaseBrowserTree::AcceptSave);
	connect(widgets->configBrowser, SIGNAL(selectionChanged()),
			this, SLOT(onSelectionChanged()));
	connect(widgets->runBrowser, SIGNAL(selectionChanged()),
			this, SLOT(onSelectionChanged()));
	connect(widgets->addButton, SIGNAL(clicked()), this, SLOT(addRunConditions()));
}

/**
 * Destructor.
 */
DbRunEditorDialog::~DbRunEditorDialog()
{
}

/**
 * Defines the model for the database browser.
 */
void DbRunEditorDialog::setDatabaseModel(Model::DatabaseModel* model)
{
	widgets->configBrowser->setDatabaseModel(model);
	widgets->runBrowser->setDatabaseModel(model);
}

/**
 * Returns the list of selected configuration versions items.
 */
std::vector<Model::AbstractDatabaseItem*> DbRunEditorDialog::getConfigSelection()
{
	return widgets->configBrowser->getSelection();
}

/**
 * Returns the list of selected run items.
 */
std::vector<Model::AbstractDatabaseItem*> DbRunEditorDialog::getRunSelection()
{
	return widgets->runBrowser->getSelection();
}

/**
 * Enables accept button when user selected a configuration and a run.
 */
void DbRunEditorDialog::onSelectionChanged()
{
    	std::vector<Model::AbstractDatabaseItem*> configSelection = getConfigSelection();
		std::vector<Model::AbstractDatabaseItem*> runSelection = getRunSelection();
    	widgets->addButton->setEnabled(not (configSelection.empty() or runSelection.empty()));
}

/**
 * Adds selected configuration versions to selected runs as conditions.
 */
void DbRunEditorDialog::addRunConditions()
{
	std::vector<Model::AbstractDatabaseItem*> configSelection = getConfigSelection();
	std::vector<Model::AbstractDatabaseItem*> runSelection = getRunSelection();
	for (std::vector<Model::AbstractDatabaseItem*>::iterator curRunItem = runSelection.begin(),
		endRunItem = runSelection.end(); curRunItem != endRunItem; ++curRunItem)
	{
		Model::RunItem* run = dynamic_cast<Model::RunItem*>(*curRunItem);
		for (std::vector<Model::AbstractDatabaseItem*>::iterator curVersionItem = configSelection.begin(),
			endVersionItem = configSelection.end(); curVersionItem != endVersionItem; ++curVersionItem)
		{
			Model::VersionItem* version = dynamic_cast<Model::VersionItem*>(*curVersionItem);
			widgets->runBrowser->getDatabaseModel()->createCondition(run, version);
		}
	}

}

}
}
