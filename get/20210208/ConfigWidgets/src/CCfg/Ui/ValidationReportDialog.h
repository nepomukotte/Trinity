#ifndef CCFG_UI_VALIDATIONREPORTDIALOG_H_INCLUDED
#define CCFG_UI_VALIDATIONREPORTDIALOG_H_INCLUDED

#include <Utils/Ui/DesignedWidget.hpp>
#include <CCfg/ValidationVisitor.h>
#include <QDialog>

// forward declarations
class Ui_ValidationReportDialog;
class QAbstractButton;

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is a dialog box implementing the configuration validation callbacks to
 * display validation errors in a list.
 *
 * @author fchateau
 */
class ValidationReportDialog
	: public Utils::Ui::DesignedWidget<Ui_ValidationReportDialog, QDialog>,
	  public AbstractValidationVisitor
{
	Q_OBJECT
public:
	// enums
	enum UserAction {Abort, Retry, Ignore, Undefined};

	ValidationReportDialog(QWidget* parent = 0);
	virtual ~ValidationReportDialog();

	bool hasErrors() const { return hasErrs; }
	UserAction getUserAction() const { return action; }

	void reset();

protected slots:
	void onButtonClicked(QAbstractButton* button);

protected:
	void onDoublet(const Attribute& attribute);
	void onOutOfRange(const Scalar& scalar);
	void onTypeMismatch(const Attribute& attr, const Attribute& other);

	bool hasErrs;
	UserAction action;

private:
	ValidationReportDialog(const ValidationReportDialog& rvalue);
	ValidationReportDialog& operator=(const ValidationReportDialog& rvalue);
};

}
}

#endif // CCFG_UI_VALIDATIONREPORTDIALOG_H_INCLUDED
