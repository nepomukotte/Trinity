#include "InteractiveUpdateController.h"
#include "UpdateControllerDialog.h"
#include "CCfg/StructureNode.h"
#include "CCfg/TypeInfo.h"
#include <QObject>

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::InteractiveUpdateController class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 * @param dialogParent the parent of the dialog, created to handle the update
 * control interactively. This parent will be responsible for freeing the dialog,
 * so it must not be null.
 */
InteractiveUpdateController::InteractiveUpdateController(QWidget* dialogParent)
	: dialog(new UpdateControllerDialog(dialogParent))
{
}

/*
 * (see inherited doc)
 */
void InteractiveUpdateController::onClassAdded(const StructureNode& newClass)
{
	QString message(QObject::tr("Class added: \"%1\" (and children classes if any)")
			.arg(newClass.path.c_str()));
	dialog->addMessage(message);
}

/*
 * (see inherited doc)
 */
void InteractiveUpdateController::onClassRemoved(const StructureNode& oldClass)
{
	QString message(QObject::tr("Class removed: \"%1\"")
			.arg(oldClass.path.c_str()));
	dialog->addMessage(message);
}

/*
 * (see inherited doc)
 */
void InteractiveUpdateController::onUnitChanged(const StructureNode& newClass, const StructureNode& oldClass)
{
	QString message(QObject::tr("Unit changed: \"%1\" switched from \"%2\" to \"%3\"")
			.arg(newClass.path.c_str())
			.arg(oldClass.unit.c_str())
			.arg(newClass.unit.c_str()));
	dialog->addMessage(message);
}

/*
 * (see inherited doc)
 */
void InteractiveUpdateController::onRangeChanged(const StructureNode& newClass, const StructureNode& oldClass)
{
	QString message(QObject::tr("Range changed: \"%1\" switched from \"%2\" to \"%3\"")
			.arg(newClass.path.c_str())
			.arg(oldClass.range.c_str())
			.arg(newClass.range.c_str()));
	dialog->addMessage(message);
}

/*
 * (see inherited doc)
 */
void InteractiveUpdateController::onTypeChanged(const StructureNode& newClass, const StructureNode& oldClass)
{
	QString message(QObject::tr("Type changed: \"%1\" switched from \"%2\" to \"%3\"")
			.arg(newClass.path.c_str())
			.arg(TypeInfo::getDisplayName(oldClass.type))
			.arg(TypeInfo::getDisplayName(newClass.type)));
	dialog->addMessage(message);
}

/*
 * (see inherited doc)
 */
bool InteractiveUpdateController::doConfirmStructureChange()
{
	int answer = dialog->exec();
	dialog->reset();
	return answer == QDialog::Accepted;
}

}
}
