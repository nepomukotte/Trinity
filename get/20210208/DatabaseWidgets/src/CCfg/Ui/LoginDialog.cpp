#include "LoginDialog.h"
#include "CCfg/Ui/LoginDialog.ui.h"

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::LoginDialog class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
LoginDialog::LoginDialog(QWidget* parent) : QDialog(parent),
	widgets(new Ui_LoginDialog)
{
	widgets->setupUi(this);
	widgets->lineEditUserName->setFocus();
}

/**
 * Destructor.
 */
LoginDialog::~LoginDialog()
{
}

/**
 * Sets the user name displayed in the dialog box.
 */
void LoginDialog::setUser(const std::string& user)
{
	widgets->lineEditUserName->setText(user.c_str());
	if(user.empty())
	{
		widgets->lineEditUserName->setFocus();
	}
	else
	{
		widgets->lineEditPassword->setFocus();
	}
}

/**
 * Retrieves the entered user name.
 */
std::string LoginDialog::getUser() const
{
	return widgets->lineEditUserName->text().toStdString();
}

/**
 * Retrieves the entered password.
 */
std::string LoginDialog::getPassword() const
{
	return widgets->lineEditPassword->text().toStdString();
}

/**
 * Clears login informations.
 */
void LoginDialog::clear()
{
	widgets->lineEditUserName->clear();
	widgets->lineEditPassword->clear();
	widgets->lineEditUserName->setFocus();
}

}
}
