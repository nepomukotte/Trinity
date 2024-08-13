#ifndef CCFG_UI_DBCONFIGCHOOSERDIALOG_H_INCLUDED
#define CCFG_UI_DBCONFIGCHOOSERDIALOG_H_INCLUDED
#include "DatabaseBrowserTree.h"
#include "CCfg/Model/AbstractDatabaseItem.h"
#include "CCfg/Model/DatabaseFilterProxyModel.h"
#include <QDialog>
#include <memory>

// forward declarations
class Ui_DbConfigChooserDialog;

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
 * This dialog box is like a QFileDialog excepts that it selects a configuration version
 * in a database.
 *
 * @author fchateau
 */
class DbConfigChooserDialog : public QDialog
{
	Q_OBJECT
public:
	DbConfigChooserDialog(QWidget* parent = 0);
	~DbConfigChooserDialog();

	void setDatabaseModel(Model::DatabaseModel* model);
	void setAcceptMode(DatabaseBrowserTree::AcceptMode newMode);
	Model::DatabaseModel* getDatabaseModel() const;
	DatabaseBrowserTree::AcceptMode getAcceptMode() const;
	void setFilterMode(DatabaseFilterProxyModel::FilterMode newMode);
	std::vector<Model::AbstractDatabaseItem*> getSelection();

protected:
	using QObject::connect;
protected slots:
	void onCurrentItemChanged(AbstractDatabaseItem* currentItem);
	void onDoubleClicked(AbstractDatabaseItem* item);
private:
	DbConfigChooserDialog(const DbConfigChooserDialog& rvalue);
	DbConfigChooserDialog& operator=(const DbConfigChooserDialog& rvalue);

	std::auto_ptr<Ui_DbConfigChooserDialog> widgets;
};

}
}

#endif // CCFG_UI_DBCONFIGCHOOSERDIALOG_H_INCLUDED
