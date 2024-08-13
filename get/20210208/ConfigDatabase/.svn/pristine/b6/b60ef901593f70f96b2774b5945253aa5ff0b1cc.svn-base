// $Id$
#include "RunCondition.h"
#include "ConfigDao.h"
#include <cassert>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::RunCondition class methods body
////////////////////////////////////////////////////////////////////////////////

const char RunCondition::TABLE_NAME[] = "CCfg_RunCondition";

/**
 * Default constructor.
 */
RunCondition::RunCondition()
{
}

/**
 * Destructor.
 */
RunCondition::~RunCondition()
{
}

/**
 * Equality operator.
 */
bool RunCondition::operator==(const RunCondition& entity)
{
	return id == entity.id;
}

/**
 * @fn bool RunCondition::operator!=(const RunCondition& entity) { return not operator==(entity); }
 * Unequality operator.
 */

/**
 * Retrieves the run associated with this run condition, from the database.
 */
RunPtr RunCondition::getRun()
{
	assert(dao != 0);
	return dao->getRun(fkRun);
}

/**
 * Retrieves the version of the configuration associated with this run condition, from the database.
 */
ConfigVersionPtr RunCondition::getConfigVersion()
{
	assert(dao != 0);
	return dao->getVersion(fkConfigVersion);
}

}
}
}
