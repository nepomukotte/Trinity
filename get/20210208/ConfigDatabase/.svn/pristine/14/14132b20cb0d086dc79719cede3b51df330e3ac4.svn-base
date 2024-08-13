// $Id$
#ifndef CCFG_IO_DB_SITE_H_INCLUDED
#define CCFG_IO_DB_SITE_H_INCLUDED

#include "ConfigEntity.h"
#include "EntitySets.h"
#include "CCfg/Attribute.h"
#include "CCfg/Io/Db/Run.h"
#include <boost/intrusive_ptr.hpp>
#include <string>
#include <vector>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// forward declarations
////////////////////////////////////////////////////////////////////////////////
struct Site;
typedef boost::intrusive_ptr<Site> SitePtr;

////////////////////////////////////////////////////////////////////////////////
// class declarations
////////////////////////////////////////////////////////////////////////////////
/**
 * This entity stands for a site with which runs and their configurations can be associated.
 * Each time a configuration is associated with a run at a new site, this site is saved into the database.
 *
 * @author psizun
 */
struct Site: public ConfigEntity
{
	static const char TABLE_NAME[];
	enum ColumnIndexes {ID = 0, NAME};

	std::string name;

	Site();
	virtual ~Site();

	UID getId() const { return id; }
	const std::string& getName() const { return name; }

	//bool operator==(const Site& entity);
	//bool operator!=(const Site& entity) { return not operator==(entity); }

	static SitePtr createSite(const std::string& siteName, ConfigDao* dao);
	RunPtr createNewRun();
	RunPtrSetPtr getRuns();

private:
	Site(const Site& rvalue);
	Site& operator=(const Site& rvalue);
};

typedef boost::multi_index_container<
	SitePtr,
	indexed_by<
		ordered_unique<tag<Id>, BOOST_MULTI_INDEX_CONST_MEM_FUN(Site, UID, getId)>,
		ordered_unique<tag<Name>, BOOST_MULTI_INDEX_CONST_MEM_FUN(Site, const std::string&, getName)>
	>
> SitePtrSet;
typedef std::auto_ptr<SitePtrSet>	SitePtrSetPtr;
typedef SitePtrSet::index<Id>::type SiteIdIdx;
typedef SitePtrSet::index<Name>::type SiteNameIdx;

}
}
}

#endif /* CCFG_IO_DB_SITE_H_INCLUDED */
