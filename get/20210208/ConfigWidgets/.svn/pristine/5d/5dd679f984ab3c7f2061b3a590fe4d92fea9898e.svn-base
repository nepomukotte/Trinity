#include "Workspace.h"
#include "ConfigEditorController.h"
#include "ReadOnlyConfigEditorController.h"
#include "../Ui/ConfigEditorTree.h"
#include "../Model/ConfigModel.h"
#include "../Model/ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include "CCfg/Io/Document.h"
#include "CCfg/Io/Url.h"
#include "../Ui/ClosableTabWidget.h"
#include "../Ui/AttributeBuilderWidget.h"
#include "../Ui/ValidationReportDialog.h"
#include "../Model/AttributeBuilderModel.h"
#include "CCfg/Object.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDialog>
#include <QVBoxLayout>
#include <QWidget>
#include <QModelIndex>
#include <QDir>
#include <QUndoStack>
#include <sstream>

namespace CCfg
{
using Model::ConfigModel;
using Model::ConfigNode;
using Ui::ConfigEditorTree;
using Ui::AttributeBuilderWidget;

namespace Control
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Control::Workspace class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param workspaceUi the tab widget that this class will take control.
 * @param parent the parent widget.
 */
Workspace::Workspace(Ui::ClosableTabWidget* workspaceUi, QWidget* parent)
	: QObject(parent),
	attributeBuilder(0), workspaceUi(0),
	controllers(), models(), currentController(0)
{
	setWorkspaceUi(workspaceUi);
}

/**
 * Destructor.
 */
Workspace::~Workspace()
{
	foreach(BaseConfigEditorController* controller, controllers)
	{
		delete controller;
	}
	controllers.clear();
	foreach(ConfigModel* model, models)
	{
		delete model;
	}
	models.clear();
}

/**
 * Defines the external attribute builder to use when a document needs to create an attribute.
 * If this function is not called by the workspace creator, the default
 * attribute builder (a dialog box) will be created when needed.
 */
void Workspace::setAttributeBuilder(Ui::AttributeBuilderWidget* builder)
{
	// if there already were an attribute builder, disconnect its signals from
	// the workspace.
	if(attributeBuilder != 0)
	{
		disconnect(attributeBuilder, 0,	this, 0);
	}
	// connect the new attribute builder.
	connect(builder,
		SIGNAL(insertClicked(Model::AttributeBuilderModel*, QString, bool)),
		this,
		SLOT(onInsertClicked(Model::AttributeBuilderModel*, QString, bool)));

	attributeBuilder = builder;
}

/**
 * Returns the attribute builder used when a document needs to create an attribute.
 * If the client have not previously defined it, we creates a default attribute
 * builder.
 */
Ui::AttributeBuilderWidget* Workspace::getAttributeBuilder()
{
	// if the are no attribute builder, we create a default one.
	if(attributeBuilder == 0)
	{
		// create a dialog box.
		QDialog* builderDialog = new QDialog(workspaceUi);
		QVBoxLayout* layout = new QVBoxLayout(builderDialog);

		// put an attribute builder widget inside it.
		Ui::AttributeBuilderWidget* builder = new Ui::AttributeBuilderWidget(builderDialog);
		layout->addWidget(builder);
		layout->setAlignment(builder, Qt::AlignTop);
		setAttributeBuilder(builder);
	}
	return attributeBuilder;
}

/**
 * Defines the workspace user interface controlled by this class.
 */
void Workspace::setWorkspaceUi(Ui::ClosableTabWidget* widget)
{
	if(workspaceUi != 0)
	{
		disconnect(workspaceUi, 0, this, 0);
	}
	connect(widget, SIGNAL(currentChanged(int)), this, SLOT(onCurrentTabChanged(int)));
	connect(widget, SIGNAL(tabClosing(int)), this, SLOT(onTabClosing(int)));
	connect(widget, SIGNAL(tabClosed(int)), this, SLOT(onTabClosed(int)));
	connect(widget, SIGNAL(canClose(int, bool*)), this, SLOT(canClose(int, bool*)));
	connect(widget, SIGNAL(canCloseAll(bool*)), this, SLOT(canCloseAll(bool*)));

	workspaceUi = widget;
}

/**
 * Returns the workspace user interface controlled by this class.
 */
Ui::ClosableTabWidget* Workspace::getWorkspaceUi()
{
	return workspaceUi;
}

/**
 * Creates a new blank configuration.
 */
BaseConfigEditorController* Workspace::newConfig()
{
	std::auto_ptr<Io::Document> doc(new Io::Document);

	ConfigModel* model = createModel(0);
	model->setDocument(doc);
	addModel(model);

	return createViewOnModel(model);
}

/**
 * Opens the files the user selects in the dialog box.
 */
void Workspace::open()
{
	openFiles(selectFilesToOpen());
}

/**
 * Saves the active model.
 * @return Returns true if the document was saved, and false if the user cancelled
 * the "save as" dialog (only occurs if the document was unnamed).
 */
bool Workspace::save()
{
	return save(getCurrentModel());
}

/**
 * Saves the active model in a new file.
 * @return Returns true if the document was saved, and false if the user cancelled
 * the "save as" dialog.
 */
bool Workspace::saveAs()
{
	return saveAs(getCurrentModel());
}

/**
 * Saves all modified configurations.
 */
void Workspace::saveAll()
{
	foreach(Model::ConfigModel* model, models)
	{
		if(not model->getUndoStack()->isClean())
		{
			save(model);
		}
	}
}

/**
 * Closes the active tab.
 * The user can interact if the document is not saved, and oppose to the closure.
 * @return Returns true if the tab was closed, false otherwise.
 */
bool Workspace::close()
{
	return removeTab(getCurrentController(), true);
}

/**
 * Closes all tabs.
 * The user can if at least one document is not saved, and oppose to the closure.
 * @return Returns true if all tabs were closed, false otherwise.
 */
bool Workspace::closeAll()
{
	bool isOk = okToCloseAll();
	if(isOk)
	{
		workspaceUi->removeAll();
	}
	return isOk;
}

/**
 * Prompts the user about what he wants to do before closing all tabs.
 * If the user chooses to save a model or all of them, the action is carried out
 * immediately.
 * @return Returns true if the user agrees to close all tab (no matter if he
 * saved them or not), false if he refuses.
 */
bool Workspace::okToCloseAll()
{
	if (1 == models.size())
	{
		return okToClose(models.first());
	}

	foreach(Model::ConfigModel* model, models)
	{
		if(not model->getUndoStack()->isClean())
		{
			QMessageBox::StandardButton r =
				QMessageBox::warning(workspaceUi, tr("Closing all.."),
	            	tr("The document '%1' has been modified.\n"
						"What do you want to do?")
						.arg(model->getDocumentName().c_str()),
		            QMessageBox::Save | QMessageBox::SaveAll | QMessageBox::Discard |
		            QMessageBox::Close | QMessageBox::Cancel,
		            QMessageBox::Cancel);

			if(r == QMessageBox::Save)
			{
				save(model);
			}
			else if(r == QMessageBox::SaveAll)
			{
				saveAll();
				break;
			}
			else if(r == QMessageBox::Close)
			{
				break;
			}
			else if(r == QMessageBox::Cancel)
			{
				return false;
			}
		}
	}
	return true;
}

/**
 * Prompts the user about what he wants to do before closing the last tab of
 * the specified model.
 * If the user chooses to save the model, the action is carried out immediately.
 * @return Returns true if the user agrees to close the model (no matter if he
 * chose to save it or not), false if he refuses.
 */
bool Workspace::okToClose(Model::ConfigModel* model)
{
	if(not model->getUndoStack()->isClean())
	{
		QMessageBox::StandardButton r =
			QMessageBox::warning(workspaceUi, tr("Tab closing"),
            	tr("The document '%1' has been modified.\n"
					"What do you want to do?")
					.arg(model->getDocumentName().c_str()),
	            QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel,
	            QMessageBox::Cancel);
        if(r == QMessageBox::Save)
        {
            return save(model);
        }
        else if(r == QMessageBox::Cancel)
        {
            return false;
        }
	}
	return true;
}


void Workspace::openFiles(const QStringList& fileList)
{
	foreach(QString fileName, fileList)
	{
		QDir file(fileName);
		openConfig(Io::Url::fromPath(file.canonicalPath().toStdString()));
	}
}

/**
 * Opens a configuration.
 * @param configUrl the URL to the configuration.
 * @param forceNewView if the configuration is already open and this parameter
 * is true, it forces to create an additionnal view on the configuration.
 */
BaseConfigEditorController* Workspace::openConfig(const Io::Url& configUrl, bool forceNewView)
{
	ConfigModel* model = searchModel(configUrl.format());
	if(model)
	{
		if(forceNewView)
		{
			return createViewOnModel(model);
		}
		else
		{
			return focusModel(model);
		}
	}
	else
	{
		return createViewForUrl(configUrl);
	}
}

/**
 * Closes a configuration view.
 * The user is not prompted and can oppose to this action.
 * @param controller the controller of the configuration view to close.
 */
void Workspace::closeConfig(BaseConfigEditorController* controller)
{
	removeTab(controller, false);
}

/**
 * Saves a model.
 * @return Returns true if the document was saved, and false if the user cancelled
 * the "save as" dialog (only occurs if the document was unnamed).
 */
bool Workspace::save(Model::ConfigModel* model)
{
	if(model->getDocumentUrl().empty())
	{
		return saveAs(model);
	}
	else
	{
		model->save();
		return true;
	}
}

/**
 * Saves a model in the file selected by the user.
 * @return Returns true if the document was saved, and false if the user cancelled
 * the "save as" dialog.
 */
bool Workspace::saveAs(Model::ConfigModel* model)
{
	QString fileName = selectFileToSave();
	if(not fileName.isEmpty())
	{
		model->saveAs(Io::Url::fromPath(fileName.toStdString()));
		return true;
	}
	return false;
}

/**
 * Put the focus on the first tab displaying the specified model.
 * @return Returns the controller of the model that has been focused.
 */
BaseConfigEditorController* Workspace::focusModel(const Model::ConfigModel* model)
{
	BaseConfigEditorController* controller = findFirstControllerForModel(model);
	assert(controller != 0);
	workspaceUi->setCurrentWidget(controller->getView());
	return controller;
}

/**
 * Creates a tab for the specified model (that may be already displayed by another tab).
 * @return Returns the controller for the new tab.
 */
BaseConfigEditorController* Workspace::createViewOnModel(Model::ConfigModel* model)
{
	BaseConfigEditorController* controller = createController(this);
	controller->setModel(model);
	controller->setView(createView(0));
	addTab(controller);
	return controller;
}

/**
 * Creates a tab for a new configuration to load.
 * @return Returns the controller for the new tab.
 */
BaseConfigEditorController* Workspace::createViewForUrl(const Io::Url& configUrl)
{
	std::auto_ptr<Io::Document> doc(new Io::Document);
	Ui::ValidationReportDialog validator(workspaceUi);

	do
	{
		// load the document
		doc->loadFromUrl(configUrl);

		// validate the config
		validator.reset();
		doc->getConfig()->accept(validator);

		if(validator.hasErrors())
		{
			validator.exec();
			if(validator.getUserAction() == Ui::ValidationReportDialog::Abort)
			{
				return 0;
			}
			else if(validator.getUserAction() == Ui::ValidationReportDialog::Ignore)
			{
				break;
			}
		}
	}
	while(validator.hasErrors() and
			validator.getUserAction() == Ui::ValidationReportDialog::Retry);

	// create the model and associate it with the document
	ConfigModel* model = createModel(0);
	model->setDocument(doc);
	addModel(model);

	// create the controller and associate it with the model and the view
	return createViewOnModel(model);
}

/**
 * Adds a new controller and creates a tab for it.
 */
void Workspace::addTab(BaseConfigEditorController* controller)
{
	ConfigModel* model = controller->getModel();
	ConfigEditorTree* view = controller->getView();

	emit tabAdded(controller);
	addController(controller);

	int tabIndex = workspaceUi->addTab(view, QString::fromStdString(model->getDocumentName()));
	workspaceUi->setTabToolTip(tabIndex, QString::fromStdString(model->getDocumentUrl()));
	workspaceUi->setCurrentIndex(tabIndex);

	view->expand(view->getConfigRoot());
}

/**
 * Removes the tab of the specified controller.
 * @param controller the controller to remove.
 * @param interactive if true, it allows user interaction before removing the
 * tab (for saving or refusing the action).
 */
bool Workspace::removeTab(BaseConfigEditorController* controller, bool interactive)
{
	return workspaceUi->removeTab(controller->getView(), interactive);
}

/**
 * Changes the ability to modify things in all configuration editors for the specified model.
 */
/*void Workspace::setModelReadOnly(Model::ConfigModel* model, bool isReadOnly)
{
	for(int i = 0; i < controllers.size(); ++i)
	{
		if(controllers[i]->getModel() == model)
		{
			ConfigEditorTree* view = controllers[i]->getView();

			BaseConfigEditorController* oldController = controllers[i];
			BaseConfigEditorController* newController =
				getFactory()->createController(this, isReadOnly);

			newController->setModel(model);
			newController->setView(view);

			if(currentController == oldController)
			{
				currentController = newController;
			}
			controllers[i] = newController;
			delete oldController;
		}
	}
}*/

/**
 * Returns the map between configuration URLs and models.
 */
const ModelsList& Workspace::getModels() const
{
	return models;
}

/**
 * Returns the list of the controllers displayed by the tabs.
 */
const ControllersList& Workspace::getControllers() const
{
	return controllers;
}

/**
 * Returns the controller of the active tab.
 */
BaseConfigEditorController* Workspace::getCurrentController()
{
	return currentController;
}

/**
 * Return the model of the active tab.
 */
Model::ConfigModel* Workspace::getCurrentModel()
{
	return currentController ? currentController->getModel() : 0;
}

/**
 * Returns the first controller (using the order of tabs) for a model.
 * @return Returns the controller if found, and 0 otherwise.
 */
BaseConfigEditorController* Workspace::findFirstControllerForModel(const Model::ConfigModel* model)
{
	foreach(BaseConfigEditorController* controller, controllers)
	{
		if(controller->getModel() == model)
		{
			return controller;
		}
	}
	return 0;
}

/**
 * Returns the total number of controllers referencing a model.
 */
unsigned Workspace::countControllersForModel(const Model::ConfigModel* model)
{
	unsigned count = 0;
	foreach(BaseConfigEditorController* controller, controllers)
	{
		if(controller->getModel() == model)
		{
			++count;
		}
	}
	return count;
}

/**
 * Searches a model matching the URL.
 * @return Returns the first model matching the URL, and 0 if no one matches it.
 */
Model::ConfigModel* Workspace::searchModel(const std::string& url)
{
	foreach(ConfigModel* model, models)
	{
		if(model->getDocumentUrl() == url)
		{
			return model;
		}
	}
	return 0;
}

/**
 * Returns the tab index of a controller.
 */
int Workspace::getControllerTabIndex(BaseConfigEditorController* controller)
{
	return workspaceUi->indexOf(controller->getView());
}

/**
 * Returns the controller of the specified view widget.
 */
BaseConfigEditorController* Workspace::getControllerForView(const QWidget* view)
{
	foreach(BaseConfigEditorController* controller, controllers)
	{
		if(controller->getView() == view)
		{
			return controller;
		}
	}
	return 0;
}

/**
 * @fn void Workspace::tabAdded(Control::BaseConfigEditorController* controller)
 * This signal is emitted when a new tab is created.
 */

/**
 * @fn void Workspace::tabAboutToBeRemoved(Control::BaseConfigEditorController* controller)
 * This signal is emitted before a tab is removed.
 */

/**
 * @fn void Workspace::tabRemoved(Control::BaseConfigEditorController* controller)
 * This signal is emitted after a tab has been removed.
 */

/**
 * @fn void Workspace::tabActivated(Control::BaseConfigEditorController* current, Control::BaseConfigEditorController* previous)
 * This signal notifies that the currently activated tab has changed.
 * It is emitted when the user selects another tab, when a tab is created
 * (because it gets activated), when a tab is removed (because another one gets
 * activated), and even when the last tab is removed (in this case "current" is
 * null).
 */

/**
 * Creates an instance of model.
 * @param parent the parent for clean-up.
 */
Model::ConfigModel* Workspace::createModel(QObject* parent)
{
	return new ConfigModel(parent);
}

/**
 * Creates an instance of view.
 * @param parent the parent for clean-up.
 */
Ui::ConfigEditorTree* Workspace::createView(QWidget* parent)
{
	return new Ui::ConfigEditorTree(parent);
}

/**
 * Creates an instance of controller.
 * @param parent the parent for clean-up.
 */
BaseConfigEditorController* Workspace::createController(QObject* parent)
{
	return new ConfigEditorController(parent);
}

/**
 * Adds a model to the list of models.
 */
void Workspace::addModel(Model::ConfigModel* model)
{
	models.append(model);
	connect(model, SIGNAL(documentUrlChanged(std::string, std::string)),
		this, SLOT(onDocumentUrlChanged(std::string, std::string)));
	emit modelAdded(model);
}

/**
 * Adds a controller to the list of controllers.
 */
void Workspace::addController(BaseConfigEditorController* controller)
{
	controllers.append(controller);
}

/**
 * Removes a model from the list of models and deletes it.
 */
void Workspace::removeModel(Model::ConfigModel* model)
{
	models.removeAll(model);
	if(models.size() == 0)
	{
		emit modelActivated(0, model);
		emit currentConfigurationChanged(0);
	}
	emit modelRemoved(model);
	delete model;
}

/**
 * Removes a controller from the list of controllers and deletes it.
 */
void Workspace::removeController(BaseConfigEditorController* controller)
{
	ConfigModel* model = controller->getModel();

	emit tabRemoved(controller);

	// removes the view
	asyncDeleteView(controller->getView());

	// removes the controller
	controllers.removeAll(controller);
	delete controller;
	if(currentController == controller)
	{
		currentController = 0;
	}

	// removes the model if necessary
	BaseConfigEditorController* otherController = findFirstControllerForModel(model);
	if(otherController == 0)
	{
		// the model is no longer referenced
		removeModel(model);
	}
}

/**
 * Does an APC (asynchronous procedure call) to delete a view.
 * The APC is necessary because when deleting the view in response to an event
 * like tabClosed, the view is still used by the functions below the signal
 * emission in the call stack.
 */
void Workspace::asyncDeleteView(Ui::ConfigEditorTree* view)
{
	QMetaObject::invokeMethod(this, "deleteView", Qt::QueuedConnection,
		Q_ARG(QWidget*, view));
}

/**
 * Displays a "file dialog" and returns the list of files to open that the user
 * has selected.
 */
QStringList Workspace::selectFilesToOpen()
{
 	return QFileDialog::getOpenFileNames(workspaceUi,
 		QString(),
 		QDir::currentPath(),
 		ConfigEditorController::FILE_EXTENSIONS_FILTER);
}

/**
 * Displays a "file dialog" and returns the list of files to save that the user
 * has selected.
 */
QString Workspace::selectFileToSave()
{
 	return QFileDialog::getSaveFileName(workspaceUi,
 		QString(),
 		QDir::currentPath(),
 		ConfigEditorController::FILE_EXTENSIONS_FILTER);
}

/**
 * Updates the target of the attribute builder.
 * @param node the new target node.
 */
void Workspace::updateAttributeBuilderTarget(Model::ConfigNode* node)
{
	if(node != 0)
	{
		getAttributeBuilder()->onTargetChanged(
			QString::fromStdString(node->getPath().toString()),
			node->isObject(), node->isRoot(), node->isInherited());
	}
	else
	{
		getAttributeBuilder()->onTargetChanged("Invalid", false, true, true);
	}
}

/**
 * Reacts to the activation of another tab.
 */
void Workspace::onCurrentTabChanged(int tabIndex)
{
	try
	{
		// When a tab has just been removed, the index cannot be used directly
		// in our controllers list, because in the QTabWidget the tab has already
		// been removed. Instead, we have to use the widget to retrieve the controller.
		BaseConfigEditorController* newController = getControllerForView(workspaceUi->widget(tabIndex));
		if(newController != 0 and currentController != newController)
		{
			ConfigModel* currentModel = currentController ? currentController->getModel() : 0;
			QItemSelectionModel* currentSelection = currentController ? currentController->getView()->selectionModel() : 0;

			ConfigModel* newModel = newController->getModel();
			QItemSelectionModel* newSelection = newController->getView()->selectionModel();

			// disconnect slots from the previous current controller.
			if(currentController != 0)
			{
				disconnect(currentSelection, 0, this, 0);
				disconnect(currentController, 0, this, 0);
			}

			// emit activation signals
			if(newModel != currentModel)
			{
				emit modelActivated(newModel, currentModel);
				emit currentConfigurationChanged(newModel->getDocument()->getConfig());
			}
			emit tabActivated(newController, currentController);

			// notify current node changed
			QModelIndex currentIndex;
			if(currentSelection != 0)
			{
				currentIndex = currentSelection->currentIndex();
			}
			QModelIndex newIndex = newSelection->currentIndex();
			ConfigNode* currentNode = currentIndex.isValid() ? currentModel->indexToNode(currentIndex) : 0;
			ConfigNode* newNode = newIndex.isValid() ? newModel->indexToNode(newIndex) : 0;
			emit currentNodeChanged(newNode, newModel, currentNode, currentModel);
			updateAttributeBuilderTarget(newNode);

			// commit controller change
			currentController = newController;

			// connect slots to the new current controller.
			connect(currentController,
				SIGNAL(attributeBuilderRequested(QPoint)),
				this,
				SLOT(onAttributeBuilderRequested(QPoint)));
			connect(currentController->getView()->selectionModel(),
				SIGNAL(currentChanged(QModelIndex, QModelIndex)),
				this,
				SLOT(onCurrentNodeChanged(QModelIndex, QModelIndex)));
		}
	}
	CATCH_AND_REPORT
}

/**
 * React to the event of a tab about to be closed.
 */
void Workspace::onTabClosing(int tabIndex)
{
	try
	{
		BaseConfigEditorController* controller = controllers[tabIndex];
		ConfigModel* model = controller->getModel();

		unsigned refCount = countControllersForModel(model);
		if(refCount == 1)
		{
			emit modelAboutToBeRemoved(model);
		}
		emit tabAboutToBeRemoved(controller);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a tab closure.
 */
void Workspace::onTabClosed(int tabIndex)
{
	try
	{
		if(controllers.size() == 1) // if it is the last tab.
		{
			assert(tabIndex == 0);
			emit tabActivated(0, controllers[tabIndex]);
		}
		removeController(controllers[tabIndex]);
	}
	CATCH_AND_REPORT
}

/**
 * Respond the TabWidget whether it is allowed to close the specified tab or not.
 */
void Workspace::canClose(int tabIndex, bool* allowed)
{
	try
	{
		*allowed = okToClose(controllers[tabIndex]->getModel());
	}
	CATCH_AND_REPORT
}

/**
 * Respond to the TabWidget whether it is allowed to close all tabs or not.
 */
void Workspace::canCloseAll(bool* allowed)
{
	try
	{
		*allowed = okToCloseAll();
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a request to display the attribute builder.
 */
void Workspace::onAttributeBuilderRequested(const QPoint& position)
{
	try
	{
		QWidget* builderParent = getAttributeBuilder()->parentWidget();
		builderParent->setFocus(Qt::PopupFocusReason);
		builderParent->show();
		if(builderParent->inherits("QDialog"))
		{
			builderParent->move(position);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a change of the current node in the current tab.
 */
void Workspace::onCurrentNodeChanged(const QModelIndex& current, const QModelIndex& previous)
{
	try
	{
		ConfigNode* currentNode = current.isValid() ? ConfigModel::indexToNode(current) : 0;
		ConfigNode* previousNode = previous.isValid() ? ConfigModel::indexToNode(previous) : 0;
		emit currentNodeChanged(currentNode,  getCurrentModel(), previousNode, getCurrentModel());
		updateAttributeBuilderTarget(currentNode);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to an insertion request from the attribute builder.
 */
void Workspace::onInsertClicked(Model::AttributeBuilderModel* model,
	const QString& path,
	bool insertAsChild)
{
	try
	{
		currentController->onInsertNew(model, path, insertAsChild);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the modification of the document URL.
 */
void Workspace::onDocumentUrlChanged(const std::string& /*newUrl*/, const std::string& /*oldUrl*/)
{
	try
	{
		ConfigModel* model = qobject_cast<ConfigModel*>(sender());
		for(int i = 0; i < controllers.size(); ++i)
		{
			if(controllers[i]->getModel() == model)
			{
				workspaceUi->setTabText(i, QString::fromStdString(model->getDocumentName()));
				workspaceUi->setTabToolTip(i, QString::fromStdString(model->getDocumentUrl()));
			}
		}
	}
	CATCH_AND_REPORT
}

/**
 * Deletes the specified view.
 * It is necessary to have this slot for asynchronous calls.
 */
void Workspace::deleteView(QWidget* view)
{
	try
	{
		delete view;
	}
	CATCH_AND_REPORT
}

}
}
