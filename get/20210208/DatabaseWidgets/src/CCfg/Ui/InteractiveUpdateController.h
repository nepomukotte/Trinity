#ifndef CCFG_UI_INTERACTIVEUPDATECONTROLLER_H_INCLUDED
#define CCFG_UI_INTERACTIVEUPDATECONTROLLER_H_INCLUDED

#include "CCfg/Io/Db/AbstractUpdateController.h"

// forward declarations
class QWidget;

namespace CCfg
{
namespace Ui
{
// forward declarations
class UpdateControllerDialog;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class spawns a dialog box to show the user the list of structure
 * modifications and let him accept them or not.
 * The dialog box lifetime is tied to the one of the parent specified in the
 * constructor.
 * If the lifetime of this class exceeds the one of the dialog, be careful that
 * no update request will be issued after the destruction of the dialog parent.
 *
 * @author fchateau
 */
class InteractiveUpdateController : public Io::Db::AbstractUpdateController
{
	UpdateControllerDialog* dialog;

public:
	InteractiveUpdateController(QWidget* dialogParent);
	virtual ~InteractiveUpdateController() {}

protected:
	void onClassAdded(const StructureNode& newClass);
	void onClassRemoved(const StructureNode& oldClass);
	void onUnitChanged(const StructureNode& newClass, const StructureNode& oldClass);
	void onRangeChanged(const StructureNode& newClass, const StructureNode& oldClass);
	void onTypeChanged(const StructureNode& newClass, const StructureNode& oldClass);
	bool doConfirmStructureChange();
};

}
}

#endif // CCFG_UI_INTERACTIVEUPDATECONTROLLER_H_INCLUDED
