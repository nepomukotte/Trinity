// $Id$
#include "Run.h"
#include "ConfigDao.h"
#include <cassert>

namespace CCfg
{
namespace Io
{
namespace Db
{
// forward declarations
struct Site;
typedef boost::intrusive_ptr<Site> SitePtr;

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::Run class methods body
////////////////////////////////////////////////////////////////////////////////

const char Run::TABLE_NAME[] = "CCfg_Run";

/**
 * Default constructor.
 */
Run::Run()
{
}

/**
 * Destructor.
 */
Run::~Run()
{
}

/**
 * Equality operator.
 */
bool Run::operator==(const Run& entity)
{
	return id == entity.id and
		number == entity.number;
}

/**
 * @fn bool Run::operator!=(const Run& entity) { return not operator==(entity); }
 * Unequality operator.
 */

/**
 * Retrieves the site associated with this run, from the database.
 */
SitePtr Run::getSite()
{
	assert(dao != 0);
	return dao->getSite(fkSite);
}

/**
 * Instantiates a new RunCondition entity for this run.
 * (It must still be inserted manually.)
 */
RunConditionPtr Run::createNewCondition(ConfigVersion& version)
{
	assert(dao != 0);
	RunConditionPtr condition(new RunCondition);
	condition->setDAO(dao);
	condition->id = dao->allocIdRangeForRunConditions(1) + 1;
	condition->fkRun = this->id;
	condition->fkConfigVersion = version.id;
	return condition;
}

/**
 * Returns the list of run conditions for this run.
 */
RunConditionPtrSetPtr Run::getConditions()
{
	assert(dao != 0);
	return dao->getConditionsOfRun(*this);
}

/**
 * Returns the list of configuration versions for this run.
 */
ConfigVersionPtrSetPtr Run::getVersions()
{
	assert(dao != 0);
	return dao->getVersionsOfRun(*this);
}

}
}
}
