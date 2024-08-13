#ifndef CCFG_UI_USERACCESSDIALOG_H_INCLUDED
#define CCFG_UI_USERACCESSDIALOG_H_INCLUDED

#include "CCfg/Io/Db/ConfigDao.h"
#include <QDialog>

// forward declarations
class Ui_UserAccessDialog;

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is a dialog box to change the user privileges on the tables of the
 * configuration system in any schema.
 *
 * @author fchateau
 */
class UserAccessDialog : public QDialog
{
	Q_OBJECT

	std::auto_ptr<Ui_UserAccessDialog> widgets;
	Io::Db::ConfigDao& dao;

public:
	UserAccessDialog(Io::Db::ConfigDao& dao, QWidget* parent = 0);
	~UserAccessDialog();

	void setSchema(const std::string& schema);
	bool isDirty() const;

public slots:
	void apply();
	void discard();

private:
	void loadUserList();
	void loadSchemaList();
	void initAccessList();

private slots:
	void onSchemaActivated(int index);
	void onUserActivated(int index);
	void onAccessActivated(int index);
	void setDirty(bool dirty);
};

}
}

#endif // CCFG_UI_USERACCESSDIALOG_H_INCLUDED
