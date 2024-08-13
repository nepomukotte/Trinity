#ifndef CCFG_IO_DB_CONFIG_H_INCLUDED
#define CCFG_IO_DB_CONFIG_H_INCLUDED

#include "ConfigStructureVersion.h"
#include "ConfigVersion.h"
#include "ConfigClass.h"
#include "EntitySets.h"
#include "ConfigEntity.h"
#include <boost/intrusive_ptr.hpp>
#include <string>
#include <vector>

namespace CCfg
{
// forward declarations
class Attribute;

namespace Io
{
namespace Db
{
// forward declarations
struct AbstractUpdateController;
struct Config;
typedef boost::intrusive_ptr<Config> ConfigPtr;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This entity stands for a configuration "file" which evolves in a new version
 * each time it is saved into the database.
 *
 * @author fchateau
 */
struct Config : public ConfigEntity
{
	static const char TABLE_NAME[];
	enum ColumnIndexes {ID = 0, NAME};

	std::string name;

	Config();
	virtual ~Config();

	UID getId() const { return id; }
	const std::string& getName() const { return name; }

	bool operator==(const Config& entity);
	bool operator!=(const Config& entity) { return not operator==(entity); }

	ConfigVersionPtr insertNewVersion(const Attribute* configRoot, AbstractUpdateController& controller);
	ConfigClassPtrSetPtr getClasses();
	ConfigStructureVersionPtrSetPtr getStructureVersions();
	ConfigVersionPtrSetPtr getVersions();
	ConfigStructureVersionPtr getLastStructureVersion();
	ConfigStructureVersionPtr getStructureVersionAt(const std::tm& date);
	ConfigVersionPtr getLastVersion();
	ConfigVersionPtr getVersionAt(const std::tm& date);
	ConfigStructureVersionPtr createNewStructure();

	static ConfigPtr createConfig(const std::string& cfgName, ConfigDao* dao);

private:
	Config(const Config& rvalue);
	Config& operator=(const Config& rvalue);
};


typedef boost::multi_index_container<
	ConfigPtr,
	indexed_by<
		ordered_unique<tag<Id>, BOOST_MULTI_INDEX_CONST_MEM_FUN(Config, UID, getId)>,
		ordered_unique<tag<Name>, BOOST_MULTI_INDEX_CONST_MEM_FUN(Config, const std::string&, getName)>
	>
> ConfigPtrSet;
typedef std::auto_ptr<ConfigPtrSet>	ConfigPtrSetPtr;
typedef ConfigPtrSet::index<Id>::type ConfigIdIdx;
typedef ConfigPtrSet::index<Name>::type ConfigNameIdx;


}
}
}

#endif // CCFG_IO_DB_CONFIG_H_INCLUDED
