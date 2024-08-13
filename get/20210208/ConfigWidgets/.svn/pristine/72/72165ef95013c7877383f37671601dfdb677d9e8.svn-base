#include "ValidationReportDialog.h"
#include "CCfg/Ui/ValidationReportDialog.ui.h"
#include "CCfg/Scalar.h"
#include <QPushButton>

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::ValidationReportDialog class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ValidationReportDialog::ValidationReportDialog(QWidget* parent)
	: Utils::Ui::DesignedWidget<Ui_ValidationReportDialog, QDialog>(parent),
	  AbstractValidationVisitor(),
	  hasErrs(false), action(Abort)
{
	connect(ui_->buttonBox, SIGNAL(clicked(QAbstractButton*)),
			this, SLOT(onButtonClicked(QAbstractButton*)));
}

/**
 * Destructor.
 */
ValidationReportDialog::~ValidationReportDialog()
{
}

/**
 * Resets the error list.
 * The object is then ready to perform a new validation.
 */
void ValidationReportDialog::reset()
{
	checkedClasses.clear();
	ui_->errorList->clear();
	hasErrs = false;
	action = Abort;
}

/**
 * Reacts to button clicks to set the action flag.
 */
void ValidationReportDialog::onButtonClicked(QAbstractButton* button)
{
	if(ui_->buttonBox->button(QDialogButtonBox::Ignore) == button)
	{
		action = Ignore;
	}
	else if(ui_->buttonBox->button(QDialogButtonBox::Retry) == button)
	{
		action = Retry;
	}
	else
	{
		action = Abort;
	}
}

/*
 * (see inherited doc)
 */
void ValidationReportDialog::onDoublet(const Attribute& attribute)
{
	hasErrs = true;
	QString message(tr("Doublet attribute!\n\"%1\" is declared more than once.")
		.arg(attribute.getPath().toString().c_str()));
	ui_->errorList->addItem(message);
}

/*
 * (see inherited doc)
 */
void ValidationReportDialog::onOutOfRange(const Scalar& /*scalar*/)
{

}

/*
 * (see inherited doc)
 */
void ValidationReportDialog::onTypeMismatch(const Attribute& attr, const Attribute& other)
{
	hasErrs = true;
	QString message(tr("Type Mismatch!\n"
		"The class \"%1\" does not have an homogenous type.\n"
		"\t\"%2\" is %3\n\t\"%4\" is %5.")
		.arg(attr.getShortClassPath().toString().c_str())
		.arg(attr.getPath().toString().c_str())
		.arg(TypeInfo::getDisplayName(attr.getType()))
		.arg(other.getPath().toString().c_str())
		.arg(TypeInfo::getDisplayName(other.getType())));
	ui_->errorList->addItem(message);
}

}
}
