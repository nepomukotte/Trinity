#ifndef CCFG_UI_DBRUNEDITORDIALOG_H_INCLUDED
#define CCFG_UI_DBRUNEDITORDIALOG_H_INCLUDED
#include "DatabaseBrowserTree.h"
#include <QDialog>
#include <memory>

// forward declarations
class Ui_DbRunEditorDialog;
class QPushButton;

namespace CCfg
{
// forward declarations
namespace Model
{
class DatabaseModel;
class AbstractDatabaseItem;
}
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This dialog box lets the user add configuration versions to a run in the database.
 *
 * @author sizun
 */
class DbRunEditorDialog : public QDialog
{
	Q_OBJECT
public:
	DbRunEditorDialog(QWidget* parent = 0);
	~DbRunEditorDialog();
	void setDatabaseModel(Model::DatabaseModel* model);

	std::vector<Model::AbstractDatabaseItem*> getConfigSelection();
	std::vector<Model::AbstractDatabaseItem*> getRunSelection();
protected:
	using QObject::connect;

protected slots:
	void onSelectionChanged();
	void addRunConditions();

private:
	DbRunEditorDialog(const DbRunEditorDialog& rvalue);
	DbRunEditorDialog& operator=(const DbRunEditorDialog& rvalue);
	std::auto_ptr<Ui_DbRunEditorDialog> widgets;
};

}
}
#endif // CCFG_UI_DBRUNEDITORDIALOG_H_INCLUDED
