#include "Workspace2.h"
#include "../Ui/DbConfigChooserDialog.h"
#include "../Ui/DbRunEditorDialog.h"
#include "../Ui/InteractiveLoginCallback.h"
#include "../Ui/InteractiveUpdateController.h"
#include "../Model/ConfigModel2.h"
#include "../Model/DatabaseModel.h"
#include "CCfg/Control/BaseConfigEditorController.h"
#include "CCfg/Ui/ClosableTabWidget.h"
#include "CCfg/Io/StorageFactory.h"
#include "CCfg/Io/FormatFactory.h"
#include "CCfg/Io/DatabaseStorageFactory.h"
#include "CCfg/Io/DatabaseFormatFactory.h"
#include "CCfg/Io/Url.h"
#include <QMessageBox>

namespace CCfg
{
namespace Control
{
using namespace Io;
using namespace Model;

////////////////////////////////////////////////////////////////////////////////
// CCfg::Control::Workspace2 class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
Workspace2::Workspace2(Ui::ClosableTabWidget* workspaceUi, QWidget* parent)
	: Workspace(workspaceUi, parent),
	loginCallback(new Ui::InteractiveLoginCallback(workspaceUi))
{
	DatabaseStorageFactory* stgFact = new DatabaseStorageFactory;
	std::auto_ptr<AbstractStorageFactory> storageFact(stgFact);
	stgFact->setLoginCallback(loginCallback);

	DatabaseFormatFactory* fmtFact = new DatabaseFormatFactory;
	std::auto_ptr<AbstractFormatFactory> formatFact(fmtFact);
	UpdateControllerPtr updateCtrl(new Ui::InteractiveUpdateController(workspaceUi));
	fmtFact->setUpdateController(updateCtrl);

	StorageFactory::instance().registerStorage(storageFact);
	FormatFactory::instance().registerFormat(formatFact);
}

/**
 * Destructor.
 */
Workspace2::~Workspace2()
{
}

/**
 * Loads a configuration or a run from a database.
 */
void Workspace2::checkout(Model::DatabaseFilterProxyModel::FilterMode mode)
{
	Model::DatabaseModel dbModel;
	dbModel.setLoginCallback(loginCallback);

	Ui::DbConfigChooserDialog dlg(getWorkspaceUi());
	dlg.setAcceptMode(Ui::DatabaseBrowserTree::AcceptLoad);
	dlg.setFilterMode(mode);
	dlg.setDatabaseModel(&dbModel);

	if (dlg.exec() == QDialog::Accepted)
	{
		std::vector<Model::AbstractDatabaseItem*> selection = dlg.getSelection();
		for(std::vector<Model::AbstractDatabaseItem*>::iterator curItem = selection.begin(),
			endItem = selection.end(); curItem != endItem; ++curItem)
		{
			Model::AbstractDatabaseItem& item = **curItem;
			std::string urlStr = item.buildUrl();
			BaseConfigEditorController* ctrl = openConfig(Io::Url(urlStr));

			// fix url stuffs
			ConfigModel2* mdl = dynamic_cast<ConfigModel2*>(ctrl->getModel());
			mdl->setDocumentUrl(std::string());
			mdl->setUrlInDatabase(urlStr);
		}
	}
}

/**
 * Loads a configuration from a database.
 */
void Workspace2::checkoutConfig()
{
	checkout(DatabaseFilterProxyModel::ConfigOnly);
}

/**
 * Saves the current configuration in the database where it comes from, or let
 * the user choose the database if it is not bound to any database.
 */
bool Workspace2::commit()
{
	return commit(getCurrentModel());
}

/**
 * Let the user choose a database configuration and saves the current
 * configuration as its new version.
 */
bool Workspace2::commitInto()
{
	return commitInto(getCurrentModel());
}

/**
 * Manage sites in the database.
 */
void Workspace2::checkoutRun()
{
	checkout(DatabaseFilterProxyModel::RunOnly);
}

/**
 * Edit runs in the database.
 */
void Workspace2::editRun()
{
	Model::DatabaseModel dbModel;
	dbModel.setLoginCallback(loginCallback);

	Ui::DbRunEditorDialog dlg(getWorkspaceUi());
	dlg.setDatabaseModel(&dbModel);

	dlg.exec();
}

/**
 * Saves the configuration in the specified model in the database where it comes
 * from, or let the user choose the database if it is not bound to any database.
 */
bool Workspace2::commit(Model::ConfigModel2* model)
{
	if(model->getUrlInDatabase().empty())
	{
		return commitInto(model);
	}
	else
	{
		model->commit();
		notifyCommitSuccess();
		return true;
	}
}

/**
 * Let the user choose a database configuration and saves the configuration in
 * the specified model as its new version.
 */
bool Workspace2::commitInto(Model::ConfigModel2* model)
{
	Model::DatabaseModel dbModel;
	dbModel.setLoginCallback(loginCallback);

	Ui::DbConfigChooserDialog dlg(getWorkspaceUi());
	dlg.setAcceptMode(Ui::DatabaseBrowserTree::AcceptSave);
	dlg.setDatabaseModel(&dbModel);

	if(dlg.exec() == QDialog::Accepted)
	{
		std::vector<Model::AbstractDatabaseItem*> selection = dlg.getSelection();
		if(selection.size() > 0)
		{
			model->commitInto(Io::Url(selection[0]->buildUrl()));
			notifyCommitSuccess();
			return true;
		}
	}
	return false;
}

/*
 * (see inherited doc)
 */
Model::ConfigModel* Workspace2::createModel(QObject* parent)
{
	return new Model::ConfigModel2(parent);
}

/**
 * Notifies the user that the configuration was committed successfully.
 */
void Workspace2::notifyCommitSuccess()
{
	QMessageBox::information(getWorkspaceUi(), tr("Database commit"),
			tr("Configuration committed successfully!"),
			QMessageBox::Close, QMessageBox::Close);
}

}
}
