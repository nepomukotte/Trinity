#include "ConfigDao.h"
#include "ConfigAttribute.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigAttribute class methods body
////////////////////////////////////////////////////////////////////////////////

const char ConfigAttribute::TABLE_NAME[] = "CCfg_ConfigAttribute";
const unsigned ConfigAttribute::INDEX_MAX_LENGTH;
const unsigned ConfigAttribute::VALUE_MAX_LENGTH;

/**
 * Default constructor.
 */
ConfigAttribute::ConfigAttribute()
{
}

/**
 * Destructor.
 */
ConfigAttribute::~ConfigAttribute()
{
}

/**
 * Equality operator.
 */
bool ConfigAttribute::operator==(const ConfigAttribute& entity)
{
	return id == entity.id and
		fkParent == entity.fkParent and
		fkVersion == entity.fkVersion and
		fkClassInfo == entity.fkClassInfo and
		index == entity.index and
		value == entity.value;
}

/**
 * @fn bool ConfigAttribute::operator!=(const ConfigAttribute& entity)
 * Unequality operator.
 */

}
}
}
