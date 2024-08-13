#ifndef CCFG_UI_WORKSPACE_H_INCLUDED
#define CCFG_UI_WORKSPACE_H_INCLUDED

#include <QObject>
#include <QList>
#include <QMap>
#include <string>

// forward declarations
class QWidget;
class QModelIndex;
class QPoint;

namespace CCfg
{
// forward declarations
class Object;
namespace Io { class Url; }
namespace Model {
	class ConfigModel;
	class AttributeBuilderModel;
	class ConfigNode;
}
namespace Ui {
	class ConfigEditorTree;
	class AttributeBuilderWidget;
	class ClosableTabWidget;
}

namespace Control
{
// forward declarations
class BaseConfigEditorController;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

typedef QList<BaseConfigEditorController*> ControllersList;
typedef QList<Model::ConfigModel*> ModelsList;

/**
 * The controller for a tabbed multidocument interface hosting configuration tree views.
 *
 * @author fchateau
 */
class Workspace : public QObject
{
	Q_OBJECT

	// @{ User interface
	Ui::AttributeBuilderWidget* attributeBuilder; ///< the attribute builder used when a document needs to create an attribute.
	Ui::ClosableTabWidget* workspaceUi;	///< the tab widget that this class controls.
	// @}

	// @{ Data Model
	ControllersList controllers;	///< the list of controllers (for each i, controllers[i] is the controller of the widget in tab i).
	ModelsList models;	///< a map associating the URL of each opened configuration to its model.
	BaseConfigEditorController* currentController;	///< the controller of the current tab.
	// @}

public:
	// inner classes
	class DefaultConfigEditorFactory;

	// @{ construction
	Workspace(Ui::ClosableTabWidget* workspaceUi, QWidget* parent = 0);
	virtual ~Workspace();
	// @}

	// @{ Compound Classes Accessors.
	void setAttributeBuilder(Ui::AttributeBuilderWidget* builder);
	Ui::AttributeBuilderWidget* getAttributeBuilder();

	void setWorkspaceUi(Ui::ClosableTabWidget* widget);
	Ui::ClosableTabWidget* getWorkspaceUi();
	// @}

	// @{ User-level Config Management
	BaseConfigEditorController* newConfig();
	void open();
	bool save();
	bool saveAs();
	void saveAll();
	bool close();
	bool closeAll();
	bool okToCloseAll();
	bool okToClose(Model::ConfigModel* model);
	// @}

	// @{ General Config Management
	void openFiles(const QStringList& fileList);
	BaseConfigEditorController* openConfig(const Io::Url& configUrl, bool forceNewView = false);
	void closeConfig(BaseConfigEditorController* controller);
	bool save(Model::ConfigModel* model);
	bool saveAs(Model::ConfigModel* model);
	// @}

	// @{ Finer Config Management
	BaseConfigEditorController* focusModel(const Model::ConfigModel* model);
	BaseConfigEditorController* createViewOnModel(Model::ConfigModel* model);
	BaseConfigEditorController* createViewForUrl(const Io::Url& configUrl);
	void addTab(BaseConfigEditorController* controller);
	bool removeTab(BaseConfigEditorController* controller, bool interactive = false);
	//void setModelReadOnly(Model::ConfigModel* model, bool isReadOnly = true);
	// @}

	// @{ Workspace Data Accessors
	const ModelsList& getModels() const;
	const ControllersList& getControllers() const;
	BaseConfigEditorController* getCurrentController();
	Model::ConfigModel* getCurrentModel();
	BaseConfigEditorController* findFirstControllerForModel(const Model::ConfigModel* model);
	unsigned countControllersForModel(const Model::ConfigModel* model);
	Model::ConfigModel* searchModel(const std::string& url);
	int getControllerTabIndex(BaseConfigEditorController* controller);
	BaseConfigEditorController* getControllerForView(const QWidget* view);
	// @}

signals:
	// @{ signals about tabs
	void tabAdded(Control::BaseConfigEditorController* controller);
	void tabAboutToBeRemoved(Control::BaseConfigEditorController* controller);
	void tabRemoved(Control::BaseConfigEditorController* controller);
	void tabActivated(Control::BaseConfigEditorController* current, Control::BaseConfigEditorController* previous);
	// @}

	// @{ signals about models
	void modelAdded(Model::ConfigModel* model);
	void modelAboutToBeRemoved(Model::ConfigModel* model);
	void modelRemoved(Model::ConfigModel* model);
	void modelActivated(Model::ConfigModel* current, Model::ConfigModel* previous);
	// @}

	// @{ Convenience signal about configuration change
	void currentConfigurationChanged(const CCfg::Object* configRoot);
	// @}

	void currentNodeChanged(Model::ConfigNode* current, Model::ConfigModel* currentModel, Model::ConfigNode* previous, Model::ConfigModel* previousModel);

protected:
	// @{ factory methods
	virtual Model::ConfigModel* createModel(QObject* parent);
	virtual Ui::ConfigEditorTree* createView(QWidget* parent);
	virtual BaseConfigEditorController* createController(QObject* parent);
	// @}
	// @{ Internal mutators
	void addModel(Model::ConfigModel* model);
	void addController(BaseConfigEditorController* controller);
	void removeModel(Model::ConfigModel* model);
	void removeController(BaseConfigEditorController* controller);
	void asyncDeleteView(Ui::ConfigEditorTree* view);
	// @}

	// @{ File Selection
	QStringList selectFilesToOpen();
	QString selectFileToSave();
	// @}

	void updateAttributeBuilderTarget(Model::ConfigNode* node);

protected slots:
	// @{ Tab related slots
	void onCurrentTabChanged(int tabIndex);
	void onTabClosing(int tabIndex);
	void onTabClosed(int tabIndex);
	void canClose(int tabIndex, bool* allowed);
	void canCloseAll(bool* allowed);
	// @}

	// @{ Configuration Editor related slots
	void onAttributeBuilderRequested(const QPoint& position);
	void onCurrentNodeChanged(const QModelIndex& current, const QModelIndex& previous);
	void onInsertClicked(Model::AttributeBuilderModel* model, const QString& path, bool insertAsChild);
	void onDocumentUrlChanged(const std::string& newUrl, const std::string& oldUrl);
	// @}

private slots:
	void deleteView(QWidget* view);
};

}
}

#endif // CCFG_UI_WORKSPACE_H_INCLUDED
