#ifndef CCFG_IO_DB_ABSTRACTUPDATECONTROLLER_H_INCLUDED
#define CCFG_IO_DB_ABSTRACTUPDATECONTROLLER_H_INCLUDED

#include "CCfg/TypeInfo.h"

namespace CCfg
{
// forward declarations
struct StructureNode;

namespace Io
{
namespace Db
{
/**
 * This interface allows to control the insertion of a configuration.
 * It receives notifications about structure changes and confirm them.
 *
 * @author fchateau
 */
struct AbstractUpdateController
{
	virtual ~AbstractUpdateController() {}
	virtual void onClassAdded(const StructureNode& newClass) = 0;
	virtual void onClassRemoved(const StructureNode& oldClass) = 0;
	virtual void onUnitChanged(const StructureNode& newClass, const StructureNode& oldClass) = 0;
	virtual void onRangeChanged(const StructureNode& newClass, const StructureNode& oldClass) = 0;
	virtual void onTypeChanged(const StructureNode& newClass, const StructureNode& oldClass) = 0;
	virtual bool doConfirmStructureChange() = 0;
};

/**
 * This is the default implementation of the update controller, it silently accept
 * any structure update.
 */
struct DefaultUpdateController : AbstractUpdateController
{
	virtual ~DefaultUpdateController() {}
	virtual void onClassAdded(const StructureNode& newClass);
	virtual void onClassRemoved(const StructureNode& oldClass);
	virtual void onUnitChanged(const StructureNode& newClass, const StructureNode& oldClass);
	virtual void onRangeChanged(const StructureNode& newClass, const StructureNode& oldClass);
	virtual void onTypeChanged(const StructureNode& newClass, const StructureNode& oldClass);
	virtual bool doConfirmStructureChange();
};

}
}
}

#endif // CCFG_IO_DB_ABSTRACTUPDATECONTROLLER_H_INCLUDED
