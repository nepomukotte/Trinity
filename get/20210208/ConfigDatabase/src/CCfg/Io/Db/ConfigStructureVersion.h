#ifndef CCFG_IO_DB_CONFIGSTRUCTUREVERSION_H_INCLUDED
#define CCFG_IO_DB_CONFIGSTRUCTUREVERSION_H_INCLUDED

#include "EntitySets.h"
#include "ConfigClass.h"
#include "ConfigClassInfo.h"
#include "ConfigVersion.h"
#include "ConfigEntity.h"
#include "CCfg/StructureNode.h"
#include <boost/intrusive_ptr.hpp>
#include <memory>
#include <ctime>

namespace CCfg
{
namespace Io
{
namespace Db
{
// forward declarations
struct Config;
typedef boost::intrusive_ptr<Config> ConfigPtr;
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This entity stands for a version of the structure of a configuration.
 * Each time an attribute is added, a unit, range or type is changed, a new
 * structure version is created.
 *
 * @author fchateau
 */
struct ConfigStructureVersion : public ConfigEntity
{
	static const char TABLE_NAME[];
	enum ColumnIndexes {ID = 0, CREATION_DATE, FK_CONFIG};

	UID fkConfig;
	std::tm creationDate;

	ConfigStructureVersion();
	virtual ~ConfigStructureVersion();

	UID getId() const { return id; }
	const std::tm& getCreationDate() const;
	time_t getCreationDateAsTimeT();

	bool operator==(const ConfigStructureVersion& entity);
	bool operator!=(const ConfigStructureVersion& entity) { return not operator==(entity); }

	ConfigClassPtrSetPtr getClasses();
	ConfigClassInfoPtrSetPtr getClassInfos();
	ConfigPtr getConfig();
	ConfigVersionPtr createNewVersion();
	std::auto_ptr<StructureNode> buildStructure(ConfigClassPtrSet& classes, ConfigClassInfoPtrSet& infos);

private:
	ConfigStructureVersion(const ConfigStructureVersion& rvalue);
	ConfigStructureVersion& operator=(const ConfigStructureVersion& rvalue);

	void initStructureNode(StructureNode& node, ConfigClass& cls, ConfigClassInfo& clsInfo);
};

typedef boost::intrusive_ptr<ConfigStructureVersion> ConfigStructureVersionPtr;
typedef boost::multi_index_container<
	ConfigStructureVersionPtr,
	indexed_by<
		ordered_unique<tag<Id>, BOOST_MULTI_INDEX_CONST_MEM_FUN(ConfigStructureVersion, UID, getId)>,
		ordered_unique<tag<Date>, BOOST_MULTI_INDEX_MEM_FUN(ConfigStructureVersion, time_t, getCreationDateAsTimeT)>
	>
> ConfigStructureVersionPtrSet;
typedef std::auto_ptr<ConfigStructureVersionPtrSet> ConfigStructureVersionPtrSetPtr;
typedef ConfigStructureVersionPtrSet::index<Id>::type ConfigStructureVersionIdIdx;
typedef ConfigStructureVersionPtrSet::index<Date>::type ConfigStructureVersionDateIdx;

}
}
}

#endif // CCFG_IO_DB_CONFIGSTRUCTUREVERSION_H_INCLUDED
