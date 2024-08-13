#ifndef CCFG_UI_LOGINDIALOG_H_INCLUDED
#define CCFG_UI_LOGINDIALOG_H_INCLUDED

#include <QDialog>
#include <memory>

// forward declarations
class Ui_LoginDialog;

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is the authentication dialog box for databases connections.
 *
 * @author fchateau
 */
class LoginDialog : public QDialog
{
	std::auto_ptr<Ui_LoginDialog> widgets;

public:
	LoginDialog(QWidget* parent = 0);
	~LoginDialog();

	void setUser(const std::string& user);
	std::string getUser() const;
	std::string getPassword() const;
	void clear();

private:
	LoginDialog(const LoginDialog& rvalue);
	LoginDialog& operator=(const LoginDialog& rvalue);
};

}
}

#endif // CCFG_UI_LOGINDIALOG_H_INCLUDED
