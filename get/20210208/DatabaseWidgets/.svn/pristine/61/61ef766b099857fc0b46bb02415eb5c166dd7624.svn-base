#include "DbConfigChooserDialog.h"
#include "CCfg/Model/AbstractDatabaseItem.h"
#include "CCfg/Model/ConfigItem.h"
#include "CCfg/Model/VersionItem.h"
#include "CCfg/Model/ConditionItem.h"
#include "CCfg/Ui/DbConfigChooserDialog.ui.h"
#include <QPushButton>

namespace CCfg
{
using namespace Model;

namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::DbConfigChooserDialog class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
DbConfigChooserDialog::DbConfigChooserDialog(QWidget* parent)
	: QDialog(parent),
	widgets(new Ui_DbConfigChooserDialog)
{
	widgets->setupUi(this);
	setFilterMode(DatabaseFilterProxyModel::ConfigOnly);
	connect(widgets->browser, SIGNAL(currentDatabaseItemChanged(AbstractDatabaseItem*)),
			this, SLOT(onCurrentItemChanged(AbstractDatabaseItem*)));
	connect(widgets->browser, SIGNAL(databaseItemDoubleClicked(AbstractDatabaseItem*)),
			this, SLOT(onDoubleClicked(AbstractDatabaseItem*)));
}

/**
 * Destructor.
 */
DbConfigChooserDialog::~DbConfigChooserDialog()
{
}

/**
 * Sets whether dialog enables loading config versions or saving a config version into a config.
 */
void DbConfigChooserDialog::setAcceptMode(DatabaseBrowserTree::AcceptMode newMode)
{
	widgets->browser->setAcceptMode(newMode);
}

/**
 * Returns whether dialog enables loading run conditions or saving a run.
 */
DatabaseBrowserTree::AcceptMode DbConfigChooserDialog::getAcceptMode() const
{
	return widgets->browser->getAcceptMode();
}

/**
 * Sets filtering mode.
 */
void DbConfigChooserDialog::setFilterMode(DatabaseFilterProxyModel::FilterMode newMode)
{
	widgets->browser->setFilterMode(newMode);
}

/**
 * Defines the model for the database browser.
 */
void DbConfigChooserDialog::setDatabaseModel(Model::DatabaseModel* model)
{
	widgets->browser->setDatabaseModel(model);
}

/**
 * Returns the DatabaseModel of the tree widget.
 */
Model::DatabaseModel* DbConfigChooserDialog::getDatabaseModel() const
{
	return widgets->browser->getDatabaseModel();
}

/**
 * Returns the list of selected items.
 */
std::vector<Model::AbstractDatabaseItem*> DbConfigChooserDialog::getSelection()
{
	return widgets->browser->getSelection();
}

/**
 * Builds current item URL and enables accept button.
 */
void DbConfigChooserDialog::onCurrentItemChanged(CCfg::Model::AbstractDatabaseItem* currentItem)
{
    if(not currentItem)
    {
    	 widgets->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);
    	 widgets->lineEditUrl->clear();
    }
    else
    {
    	if (getAcceptMode() == DatabaseBrowserTree::AcceptLoad)
    	{
    		VersionItem* version = dynamic_cast<VersionItem*>(currentItem);
    		ConditionItem* condition = dynamic_cast<ConditionItem*>(currentItem);
			widgets->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(version != 0 or condition != 0);
    	}
    	else
    	{
    		ConfigItem* config = dynamic_cast<ConfigItem*>(currentItem);
    		widgets->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(config != 0);
    	}
    	widgets->lineEditUrl->setText(currentItem->buildUrl().c_str());
    }
}

/**
 *  * Accepts the dialog box when the user double clicks on configs or versions.
 *   */
void DbConfigChooserDialog::onDoubleClicked(AbstractDatabaseItem*)
{
	accept();
}

}
}
