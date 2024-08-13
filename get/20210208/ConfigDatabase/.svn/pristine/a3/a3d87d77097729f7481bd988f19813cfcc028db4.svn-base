// $Id$
#include "ConfigDao.h"
#include <cassert>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::Site class methods body
////////////////////////////////////////////////////////////////////////////////

const char Site::TABLE_NAME[] = "CCfg_Site";

/**
 * Default constructor.
 */
Site::Site()
{
}

/**
 * Destructor.
 */
Site::~Site()
{
}

/**
 * Equality operator.
 */
/**bool Site::operator==(const Config& entity)
{
	return id == entity.id and
		name == entity.name;
}**/

/**
 * @fn bool Site::operator!=(const Site& entity) { return not operator==(entity); }
 * Unequality operator.
 */

/**
 * Creates a new Site entity (it must still be inserted manually).
 */
SitePtr Site::createSite(const std::string& siteName, ConfigDao* dao)
{
	assert(dao != 0);
	SitePtr site(new Site);
	site->setDAO(dao);
	site->id = dao->allocIdRangeForSites(1) + 1;
	site->name = siteName;
	return site;
}

/**
 * Instantiates a new Run entity for this site.
 * (It must still be inserted manually.)
 * Run number will only bet set at insertion.
 */
RunPtr Site::createNewRun()
{
	assert(dao != 0);
	RunPtr run(new Run);
	run->setDAO(dao);
	run->id = dao->allocIdRangeForRuns(1) + 1;
	run->fkSite = this->id;
	return run;
}

/**
 * Returns the list of runs for this site.
 */
RunPtrSetPtr Site::getRuns()
{
	assert(dao != 0);
	return dao->getRunsOfSite(*this);
}

}
}
}
