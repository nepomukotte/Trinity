#ifndef CCFG_IO_DB_CONFIGDATABASE_H_INCLUDED
#define CCFG_IO_DB_CONFIGDATABASE_H_INCLUDED

#include "EntitySets.h"
#include "Config.h"
#include <boost/intrusive_ptr.hpp>
#include <string>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This pseudo entity stands for a "database" (also known as "schema") inside a
 * database server.
 *
 * @author fchateau
 */
struct ConfigDatabase : public ConfigEntity
{
	std::string name;

	ConfigDatabase();
	virtual ~ConfigDatabase();

	UID getId() const { return id; }
	const std::string& getName() const;

	ConfigPtrSetPtr getConfigs();
	ConfigPtr searchConfig(const std::string& cfgName);
	bool checkHasAllTables();

private:
	ConfigDatabase(const ConfigDatabase& rvalue);
	ConfigDatabase& operator=(const ConfigDatabase& rvalue);
};

typedef boost::intrusive_ptr<ConfigDatabase> ConfigDatabasePtr;
typedef boost::multi_index_container<
	ConfigDatabasePtr,
	indexed_by<
		ordered_unique<tag<Name>, BOOST_MULTI_INDEX_CONST_MEM_FUN(ConfigDatabase, const std::string&, getName)>
	>
> ConfigDatabasePtrSet;
typedef std::auto_ptr<ConfigDatabasePtrSet>	ConfigDatabasePtrSetPtr;
typedef ConfigDatabasePtrSet::index<Name>::type ConfigDatabaseNameIdx;

}
}
}

#endif // CCFG_IO_DB_CONFIGDATABASE_H_INCLUDED
