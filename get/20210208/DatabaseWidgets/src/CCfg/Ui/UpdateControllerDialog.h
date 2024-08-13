#ifndef CCFG_UI_UPDATECONTROLLERDIALOG_H_INCLUDED
#define CCFG_UI_UPDATECONTROLLERDIALOG_H_INCLUDED

#include <QDialog>
#include <memory>

// forward declarations
class Ui_UpdateControllerDialog;

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class is a dialog box that shows messages describing structural changes
 * before a database insertion, and let the user accept or refuse them.
 *
 * @author fchateau
 */
class UpdateControllerDialog : public QDialog
{
public:
	UpdateControllerDialog(QWidget* parent = 0);
	virtual ~UpdateControllerDialog();

	void reset();
	void addMessage(const QString& msg);

protected:
	std::auto_ptr<Ui_UpdateControllerDialog> widgets;
};

}
}

#endif // CCFG_UI_UPDATECONTROLLERDIALOG_H_INCLUDED
