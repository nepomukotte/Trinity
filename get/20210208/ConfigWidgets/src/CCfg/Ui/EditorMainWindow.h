#ifndef CCFG_UI_EDITORMAINWINDOW_H_INCLUDED
#define CCFG_UI_EDITORMAINWINDOW_H_INCLUDED

#include <Utils/Ui/DesignedWidget.hpp>
#include <QMainWindow>
#include <QUndoGroup>

// forward declarations
class QDockWidget;
class Ui_MainWindow;

namespace CCfg
{
// forward declarations
namespace Model {
	class ConfigNode;
	class ConfigModel;
}
namespace Control {
	class Workspace;
	class BaseConfigEditorController;
}
namespace Ui
{

// forward declarations
class AttributeBuilderWidget;
class AttributeFinderDialog;
class ClosableTabWidget;


///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class implements the configuration editor main window.
 *
 * @author fchateau
 */
class EditorMainWindow : public Utils::Ui::DesignedWidget<Ui_MainWindow, QMainWindow>
{
	Q_OBJECT

	// @{ Workspace related objects
	Control::Workspace* workspaceCtrl;
	ClosableTabWidget* workspaceUi;
	// @}
protected:
	// @{ External widgets
	AttributeBuilderWidget* attributeBuilder;
	AttributeFinderDialog* attributeFinder;
	// @}

	/// @{ Undo
	QUndoGroup undoGroup;
	QAction* actionUndo;
    QAction* actionRedo;
    /// @}

public:
	EditorMainWindow(QWidget* parentWidget = 0);
	virtual ~EditorMainWindow();

	void setAttributeBuilder(AttributeBuilderWidget* attributeBuilder);

	bool isWorkspaceClean() const;
	Control::Workspace* getWorkspace();
	ClosableTabWidget* getWorkspaceUi();
	AttributeBuilderWidget* getAttributeBuilder();
	AttributeFinderDialog* getAttributeFinder();

	// @{ ui items accessors
	QMenuBar* getMenubar() const;
	QStatusBar* getStatusbar() const;
	QToolBar* getToolBar() const;
	QMenu* getMenuFile() const;
	QMenu* getMenuEdit() const;
	QMenu* getMenuView() const;
	QMenu* getMenuHelp() const;
	QMenu* getMenuColumns() const;
	QAction* getActionNew() const;
	QAction* getActionOpen() const;
	QAction* getActionSave() const;
    QAction* getActionSaveAs() const;
    QAction* getActionSaveAll() const;
    QAction* getActionQuit() const;
    QAction* getActionFind() const;
    QAction* getActionShowName() const;
    QAction* getActionShowValue() const;
    QAction* getActionShowUnit() const;
    QAction* getActionShowRange() const;
    QAction* getActionShowInclude() const;
    QAction* getActionShowBacktrackLevel() const;
    QAction* getActionAbout() const;
    // @}

    virtual void closeEvent(QCloseEvent *event);

public slots:
	// @{ Slots for File menu actions.
	void newConfig();
	void openFile();
	void openFiles(const QStringList& fileList);
	bool save();
	bool saveAs();
	void saveAll();
	// @}
	/// @name
	/// @{ Slots for Edit menu actions.
	void findAttribute();
	void findNextAttribute();
	/// @}
	// @{ Slots for Help menu actions.
	void about();
	// @}
	// @{ Slots for Window menu actions
	void showColumn(int index, bool show = true);
	void showValueColumn(bool show = true);
	void showUnitColumn(bool show = true);
	void showRangeColumn(bool show = true);
	void showBacktrackLevelColumn(bool show = true);
	void showIncludeUrlColumn(bool show = true);
	// @}

protected:
	virtual Control::Workspace* createWorkspace(ClosableTabWidget* ui);

	void saveWindowShape();
	void restoreWindowShape();

protected slots:
	// @{ Slots for model events
	virtual void onModelAdded(Model::ConfigModel* model);
	virtual void onModelRemoved(Model::ConfigModel* model);
	virtual void onModelActivated(Model::ConfigModel* current, Model::ConfigModel* previous);
	// @}

	virtual void onTabActivated(Control::BaseConfigEditorController* current, Control::BaseConfigEditorController* previous);
	virtual void onCleanChanged(bool isClean);
	virtual void onCurrentNodeChanged(Model::ConfigNode* current,
			Model::ConfigModel* currentModel,
			Model::ConfigNode* previous, Model::ConfigModel* previousModel);

	virtual void onSearchAttribute(const QString & searchText, bool caseSens, bool exactMatch, int column);

private:
	// @{ Initialization
	void connectActions();
	void initAttributeBuilder();
	void initAttributeFinder();
	void initWorkspace();
	void setupUndoFramework();
	// @}
};

}
}

#endif // CCFG_UI_EDITORMAINWINDOW_H_INCLUDED
