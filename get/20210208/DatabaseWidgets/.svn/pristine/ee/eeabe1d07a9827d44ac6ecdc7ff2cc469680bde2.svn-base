#ifndef CCFG_UI_EDITORMAINWINDOW2_H_INCLUDED
#define CCFG_UI_EDITORMAINWINDOW2_H_INCLUDED

#include "CCfg/Ui/EditorMainWindow.h"
#include "../Control/Workspace2.h"

// forward declarations
class QAction;
class QMenu;

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * The main window of the configuration editor, including database features.
 *
 * @author fchateau
 */
class EditorMainWindow2 : public EditorMainWindow
{
	Q_OBJECT
public:
	EditorMainWindow2(QWidget* parentWidget = 0);
	virtual ~EditorMainWindow2();

	Control::Workspace2* getWorkspace();

public slots:
	void checkout();
	void commit();
	void commitInto();
	void checkoutRun();
	void editRun();

protected:
	virtual Control::Workspace* createWorkspace(ClosableTabWidget* ui);

	QAction* actionCheckout;
	QAction* actionCommit;
	QAction* actionCommitInto;

	QMenu* menuRun;
	QAction* actionRunCheckout;
	QAction* actionRunEdit;

protected slots:
	void onTabActivated(Control::BaseConfigEditorController* current, Control::BaseConfigEditorController* previous);

private:
	// @{ Initialization
	void createMenuRun();
	// @}
	EditorMainWindow2(const EditorMainWindow2& rvalue);
	EditorMainWindow2& operator=(const EditorMainWindow2& rvalue);
};

}
}

#endif // CCFG_UI_EDITORMAINWINDOW2_H_INCLUDED
