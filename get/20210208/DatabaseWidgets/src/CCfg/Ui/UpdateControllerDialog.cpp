#include "UpdateControllerDialog.h"
#include "CCfg/Ui/UpdateControllerDialog.ui.h"

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::UpdateControllerDialog class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
UpdateControllerDialog::UpdateControllerDialog(QWidget* parent)
	: QDialog(parent),
	widgets(new Ui_UpdateControllerDialog)
{
	widgets->setupUi(this);
	widgets->buttonBox->setFocus();
}

/**
 * Destructor.
 */
UpdateControllerDialog::~UpdateControllerDialog()
{
}

/**
 * Resets the list in the dialog box.
 */
void UpdateControllerDialog::reset()
{
	widgets->list->clear();
}

/**
 * Appends a message to the list.
 */
void UpdateControllerDialog::addMessage(const QString& msg)
{
	widgets->list->addItem(msg);
}

}
}
