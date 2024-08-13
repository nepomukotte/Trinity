// $Id$
#ifndef CCFG_IO_DB_RUN_H_INCLUDED
#define CCFG_IO_DB_RUN_H_INCLUDED

#include "ConfigEntity.h"
#include "EntitySets.h"
#include "CCfg/Attribute.h"
#include "RunCondition.h"
#include "ConfigVersion.h"
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
struct Run;
typedef boost::intrusive_ptr<Run> RunPtr;
struct Site;
typedef boost::intrusive_ptr<Site> SitePtr;

////////////////////////////////////////////////////////////////////////////////
// class declarations
////////////////////////////////////////////////////////////////////////////////
/**
 * This entity stands for a run associated to a site. this run is characterized by a set of configurations.
 * Each time a set of configurations is associated to a site, a new run is saved into the database.
 *
 * @see Site
 * @author psizun
 */
struct Run: public ConfigEntity
{
	static const char TABLE_NAME[];
	enum ColumnIndexes {ID = 0, FK_SITE, NUMBER};

	UID fkSite;
	UID number;

	Run();
	virtual ~Run();

	UID getId() const { return id; }
	const UID& getNumber() const { return number; }

	bool operator==(const Run& entity);
	bool operator!=(const Run& entity) { return not operator==(entity); }

	SitePtr getSite();
	RunConditionPtr createNewCondition(ConfigVersion& version);
	RunConditionPtrSetPtr getConditions();
	ConfigVersionPtrSetPtr getVersions();

private:
	Run(const Run& rvalue);
	Run& operator=(const Run& rvalue);
};

typedef boost::multi_index_container<
	RunPtr,
	indexed_by<
		ordered_unique<tag<Id>, BOOST_MULTI_INDEX_CONST_MEM_FUN(Run, UID, getId)>
	>
> RunPtrSet;
typedef std::auto_ptr<RunPtrSet>	RunPtrSetPtr;
typedef RunPtrSet::index<Id>::type RunIdIdx;

}
}
}

#endif /* CCFG_IO_DB_RUN_H_INCLUDED */
