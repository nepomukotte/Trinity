#include "ConfigClassInfo.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigClassInfo class methods body
////////////////////////////////////////////////////////////////////////////////

const char ConfigClassInfo::TABLE_NAME[] = "CCfg_ConfigClassInfo";
const unsigned ConfigClassInfo::UNIT_MAX_LENGTH;
const unsigned ConfigClassInfo::RANGE_MAX_LENGTH;

/**
 * Default constructor.
 */
ConfigClassInfo::ConfigClassInfo()
{
}

/**
 * Destructor.
 */
ConfigClassInfo::~ConfigClassInfo()
{
}

/**
 * Returns the foreign key to the configuration class.
 */
UID ConfigClassInfo::getFkClass() const
{
	return fkClass;
}

/**
 * Equality operator.
 */
bool ConfigClassInfo::operator==(const ConfigClassInfo& entity)
{
	return id == entity.id and
		fkClass == entity.fkClass and
		fkStructure == entity.fkStructure and
		type == entity.type and
		unit == entity.unit and
		range == entity.range;
}

/**
 * @fn bool ConfigClassInfo::operator!=(const ConfigClassInfo& entity)
 * Unequality operator.
 */

}
}
}
