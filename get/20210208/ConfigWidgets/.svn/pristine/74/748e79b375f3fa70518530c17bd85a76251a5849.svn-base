#ifndef CCFG_CONTROL_CONFIGEDITORCONTROLLER_H_INCLUDED
#define CCFG_CONTROL_CONFIGEDITORCONTROLLER_H_INCLUDED

#include "BaseConfigEditorController.h"
#include <vector>
#include <memory>

// forward declarations
class QPoint;
class QAction;
class QMenu;
class QString;
class QModelIndex;
class QItemSelection;

namespace CCfg
{
// forward declarations
namespace Model { class ConfigNode; }

namespace Control
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class provides standard user interactions facilities on a configuration
 * editor tree.
 * It manage interactions for edition delegates and contextual menus.
 *
 * @author fchateau
 */
class ConfigEditorController : public BaseConfigEditorController
{
	Q_OBJECT
protected:
	std::auto_ptr<QMenu> menuRootNode;
	std::auto_ptr<QMenu> menuObjectNode;
	std::auto_ptr<QMenu> menuScalarNode;
	std::auto_ptr<QMenu> menuMultiselection;

	QMenu* menuInheritance;
	QMenu* menuOverrides;
	QMenu* menuOverridenBy;
	QMenu* menuInherits;
	QMenu* menuInheritedBy;
	QMenu* menuDefinedIn;
	QMenu* menuImport;
	QMenu* menuExport;
	QMenu* menuInclude;

	QAction* actionNew;
	QAction* actionOverride;
	QAction* actionRename;
	QAction* actionDelete;
	QAction* actionCopy;
	QAction* actionCut;
	QAction* actionPaste;
	QAction* actionSearch;
	QAction* actionExpand;
	QAction* actionCollapse;
	QAction* actionImportFile;
	QAction* actionImportDB;
	QAction* actionExportFile;
	QAction* actionExportDB;
	QAction* actionIncludeFile;
	QAction* actionIncludeDB;

	enum ActionId {ActionNew, ActionOverride, ActionRename, ActionDelete,
		ActionCopy, ActionCut, ActionPaste,
		ActionSearch, ActionExpand, ActionCollapse,
		ActionImport, ActionExport, ActionInclude};

public:

	// static constants
	static const char FILE_EXTENSIONS_FILTER[];

	// @{ construction
	ConfigEditorController(QObject* parent = 0);
	virtual ~ConfigEditorController();
	// @}

	void setView(Ui::ConfigEditorTree* newView);

signals:
	void attributeBuilderRequested(const QPoint& position);

public slots:
	void onInsertNew(Model::AttributeBuilderModel* builder, const QString& path, bool insertAsChild);

protected slots:
	void onContextMenu(const QPoint& pos);
	void onInheritanceMenuAboutToShow();
	void onNodeRefClicked();

	void onNew();
	void onOverride();
	void onRename();
	void onDelete();
	void onCopy();
	void onCut();
	void onPaste();
	void onSearch();
	void onExpand();
	void onCollapse();
	void onImportFile();
	void onImportDB();
	void onExportFile();
	void onExportDB();
	void onIncludeFile();
	void onIncludeDB();
	void onSave();
	void onCommit();

protected:
	void initActions(QWidget* parent);
	void initMenus();

	void fillMenuOverrides(Model::ConfigNode* node);
	void fillMenuOverridenBy(Model::ConfigNode* node);
	void fillMenuInherits(Model::ConfigNode* node);
	void fillMenuInheritedBy(Model::ConfigNode* node);
	void fillMenuDefinedIn(Model::ConfigNode* node);

	void enableMenuItems();

	bool isActionAuthorized(ActionId actionId, const QItemSelection& selection) const;
	bool isActionAuthorized(ActionId actionId, const QModelIndex& nodeIndex) const;

	void notifyUnauthorizedAction() const;

	QString selectInputFile();
	QString selectOutputFile();

private:
	void createRootContextMenu(QMenu& menu);
	void createMultiselectionContextMenu(QMenu& menu);
	void createAttributeContextMenu(QMenu& menu, bool isObject);
	void fillInheritanceMenu(QMenu& menu, const std::vector<Model::ConfigNode*>& list);

	bool isMultiselection(const QItemSelection& selection) const;

	// non-copyable
	ConfigEditorController(const ConfigEditorController&);
	ConfigEditorController& operator=(const ConfigEditorController&);
};

}
}

#endif // CCFG_CONTROL_CONFIGEDITORCONTROLLER_H_INCLUDED
