#include "ConfigDao.h"
#include "ConfigBatch.h"
#include "ConfigAttributeBatch.h"
#include "ConfigClassBatch.h"
#include "ConfigClassInfoBatch.h"
#include "ConfigStructureVersionBatch.h"
#include "ConfigVersionBatch.h"
#include "SiteBatch.h"
#include "RunBatch.h"
#include "RunConditionBatch.h"
#include "ConfigDatabase.h"
#include "Transaction.h"
#include "TableBuilder.h"
#include <soci.h>
#include <sstream>

namespace CCfg
{
namespace Io
{
namespace Db
{
using namespace SOCI;
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigDao class methods body
////////////////////////////////////////////////////////////////////////////////

const char ConfigDao::SEQUENCE_TABLE_NAME[] = "CCfg_ConfigSequence";
/**
 * This is the list of the names of all tables involved with configuration.
 */
const char* ConfigDao::TABLES[] =
{
	Config::TABLE_NAME,
	ConfigAttribute::TABLE_NAME,
	ConfigClass::TABLE_NAME,
	ConfigClassInfo::TABLE_NAME,
	ConfigStructureVersion::TABLE_NAME,
	ConfigVersion::TABLE_NAME,
	Site::TABLE_NAME,
	Run::TABLE_NAME,
	RunCondition::TABLE_NAME,
	SEQUENCE_TABLE_NAME
};

/**
 * This array defines the privileges associated to each couple (table, role)
 */
const unsigned ConfigDao::TABLE_PRIVILEGES[][ROLE_COUNT] =
{
	{Select | Insert | Update, Select, Select}, // privileges on Config
	{Select | Insert, Select | Insert, Select}, // on ConfigAttribute
	{Select | Insert, Select, Select},			// on ConfigClass
	{Select | Insert, Select, Select},			// on ConfigClassInfo
	{Select | Insert, Select, Select},			// on ConfigStructureVersion
	{Select | Insert | Update, Select | Insert, Select},	// on ConfigVersion
	{Select | Insert | Update, Select, Select}, // on Site
	{Select | Insert | Update, Select, Select}, // on Run
	{Select | Insert | Update, Select, Select}, // on RunCondition
	{Select | Update, Select | Update, Select},	// on ConfigSequence
};




/**
 * Default constructor.
 */
ConfigDao::ConfigDao() : session()
{
}

/**
 * Destructor.
 */
ConfigDao::~ConfigDao()
{
}

/**
 * @fn void ConfigDao::connect(const std::string& connectionString)
 * Connects to a database.
 * All details are specified in the connectionString and are database specific.
 * See SOCI connection strings for reference.
 */

/**
 * @fn bool ConfigDao::isConnected() const
 * Returns whether this DAO is connected to a database or not.
 */

/**
 * @fn ConfigDatabasePtr ConfigDao::getCurrentDatabase()
 * Return an entity representing the current database.
 */

/**
 * @fn ConfigDatabasePtr ConfigDao::searchDatabase(const std::string& dbName)
 * Returns the database entity bearing the specified name, if it exists.
 */

/**
 * @fn ConfigDatabasePtrSetPtr ConfigDao::getDatabases()
 * Returns the list of databases in a server.
 */

/**
 * @fn std::vector<std::string> ConfigDao::getUserNames()
 * Returns the list of users on the server.
 */

/**
 * @fn void ConfigDao::createDatabaseAndUser(const std::string& name, const std::string& pwd)
 * Creates a database and a user with the specified name and password.
 * A successful execution of this method will require to be logged as DBA.
 * In Oracle schemas (=database) and users are the same object, so to provide
 * compatibility, we have to create both.
 * If this is undesirable, you might consider creating databases and users
 * outside the framework, and then setup tables and privileges inside.
 */

/**
 * @fn void ConfigDao::removeDatabaseAndUser(const std::string& name)
 * Removes the database and the user with the specified name.
 * A successful execution of this method will require to be logged as DBA.
 */

/**
 * @fn void ConfigDao::changeUserPassword(const std::string& userName, const std::string& newPassword)
 * Changes the password of the specified user.
 */

/**
 * @fn int ConfigDao::getPrivileges(const std::string& userName, const std::string& tableName)
 * Returns the privileges granted to a user, on a table.
 */

/**
 * Grants a user the privileges associated to a role on a table.
 */
void ConfigDao::grantPrivileges(const std::string& userName,
		const std::string& tableName, int privileges)
{
	if(privileges != 0)
	{
		session->once <<
			"GRANT " << formatPrivileges(privileges) << " "
			"ON " << dbName << '.' << tableName << " "
			"TO " << userName;
	}
}

/**
 * Takes back from a user some privileges for each table part of the
 * configuration system in a database.
 * @param userName the user whose privileges will be revoked.
 * @param privileges a combination of PrivilegeFlags.
 */
void ConfigDao::revokePrivileges(const std::string& userName, int privileges)
{
	Transaction transac(*session);
	for(size_t i = 0; i < TableCount+1; ++i)
	{
		revokePrivileges(userName, TABLES[i], privileges);
	}
	transac.commit();
}

/**
 * Takes back from a user the privileges associated with a role on a table.
 */
void ConfigDao::revokePrivileges(const std::string& userName,
		const std::string& tableName, int privileges)
{
	if(privileges != 0)
	{
		session->once <<
			"REVOKE " << formatPrivileges(privileges) << " "
			"ON " << dbName << '.' << tableName << " "
			"FROM " << userName;
	}
}

/**
 * @fn void ConfigDao::revokePrivileges(const std::string& userName)
 * Takes back all privileges granted to a user on the configuration tables of a database.
 */

/**
 * Return the access level for the specified user.
 */
ConfigUserRole ConfigDao::getRole(const std::string& userName)
{
	ConfigUserRole access = Expert;
	for(size_t i = 0; i < TableCount+1; ++i)
	{
		access = std::max(access, getRole(userName, i));
	}
	return access;
}

/**
 * Return the access level for the specified user and table.
 */
ConfigUserRole ConfigDao::getRole(const std::string& userName, int tableIndex)
{
	unsigned access = getPrivileges(userName, TABLES[tableIndex]);
	if((TABLE_PRIVILEGES[tableIndex][Expert] & access) == TABLE_PRIVILEGES[tableIndex][Expert])
	{
		return Expert;
	}
	else if((TABLE_PRIVILEGES[tableIndex][Shifter] & access) == TABLE_PRIVILEGES[tableIndex][Shifter])
	{
		return Shifter;
	}
	else if((TABLE_PRIVILEGES[tableIndex][Analyst] & access) == TABLE_PRIVILEGES[tableIndex][Analyst])
	{
		return Analyst;
	}
	else
	{
		return NoAccess;
	}
}

/**
 * Sets the exact privileges of a user on the tables of the configuration system.
 */
void ConfigDao::setRole(const std::string& userName, ConfigUserRole role)
{
	Transaction transac(*session);
	for(size_t i = 0; i < TableCount+1; ++i)
	{
		revokePrivileges(userName, TABLES[i], getPrivileges(userName, TABLES[i]));
		if(role != NoAccess)
		{
			grantPrivileges(userName, TABLES[i], TABLE_PRIVILEGES[i][role]);
		}
	}
	transac.commit();
}

/**
 * Grants a user the privileges associated to a role for each table part of the
 * configuration system.
 * The privileges may not be the same between different tables.
 */
void ConfigDao::grantRole(const std::string& userName, ConfigUserRole role)
{
	if(role != NoAccess)
	{
		Transaction transac(*session);
		for(size_t i = 0; i < TableCount+1; ++i)
		{
			grantPrivileges(userName, TABLES[i], TABLE_PRIVILEGES[i][role]);
		}
		transac.commit();
	}
}

/**
 * @fn bool ConfigDao::isSchemaAdmin(const std::string& userName)
 * Returns whether the specified user can administrate the current schema.
 */

/**
 * @fn bool ConfigDao::checkTableExist(const char* tableName)
 * Returns whether a table with the specified name exists.
 */

/**
 * Returns whether the current database has all the required tables for the
 * configuration system.
 * Note that it does not check that the tables have the required column.
 */
bool ConfigDao::checkHasAllTables()
{
	return checkTableExist(SEQUENCE_TABLE_NAME) and
		checkTableExist(RunCondition::TABLE_NAME) and
		checkTableExist(Run::TABLE_NAME) and
		checkTableExist(Site::TABLE_NAME) and
		checkTableExist(ConfigAttribute::TABLE_NAME) and
		checkTableExist(ConfigClassInfo::TABLE_NAME) and
		checkTableExist(ConfigClass::TABLE_NAME) and
		checkTableExist(ConfigVersion::TABLE_NAME) and
		checkTableExist(ConfigStructureVersion::TABLE_NAME) and
		checkTableExist(Config::TABLE_NAME);
}

/**
 * Destroys all tables of the configuration system.
 * @param nothrow do not throw exceptions in case of failure (mostly if the tables
 * don't exist).
 */
void ConfigDao::dropAllTables(bool nothrow)
{
	Transaction transac(*session);
	dropTable(SEQUENCE_TABLE_NAME, nothrow);
	dropTable(RunCondition::TABLE_NAME, nothrow);
	dropTable(Run::TABLE_NAME, nothrow);
	dropTable(Site::TABLE_NAME, nothrow);
	dropTable(ConfigAttribute::TABLE_NAME, nothrow);
	dropTable(ConfigClassInfo::TABLE_NAME, nothrow);
	dropTable(ConfigClass::TABLE_NAME, nothrow);
	dropTable(ConfigVersion::TABLE_NAME, nothrow);
	dropTable(ConfigStructureVersion::TABLE_NAME, nothrow);
	dropTable(Config::TABLE_NAME, nothrow);
	transac.commit();
}

/**
 * Destroys the specified table.
 * @param tableName the name of the table to drop
 * @param nothrow do not throw exceptions in case of failure (mostly if the table
 * doesn't exist).
 */
void ConfigDao::dropTable(const char* tableName, bool nothrow)
{
	try
	{
		session->once << "DROP TABLE " << dbName << '.' << tableName;
	}
	catch(SOCI::SOCIError& e)
	{
		if(not nothrow)
		{
			throw;
		}
	}
}

/**
 * Re-creates all tables.
 */
void ConfigDao::replaceAllTables()
{
	dropAllTables(true);
	createAllTables();
}

/**
 * Creates all tables for the configuration system.
 */
void ConfigDao::createAllTables()
{
	Transaction transac(*session);
	createConfigTable();
	createConfigStructureVersionTable();
	createConfigVersionTable();
	createConfigClassTable();
	createConfigClassInfoTable();
	createConfigAttributeTable();
	createSiteTable();
	createRunTable();
	createRunConditionTable();
	createConfigSequenceTable();
	transac.commit();
}

/**
 * Creates the table of sequences for manual primary key generation.
 * We need to control primary key generation, because we need to insert batches
 * of entities. As these entities have foreign key relationships and even self
 * referencing foreign keys, we have to know what the id will be by advance.
 */
void ConfigDao::createConfigSequenceTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(SEQUENCE_TABLE_NAME);
	builder->addColumn(Config::TABLE_NAME, Integer, Unsigned, NotNull);
	builder->addColumn(ConfigAttribute::TABLE_NAME, Integer, Unsigned, NotNull);
	builder->addColumn(ConfigClass::TABLE_NAME, Integer, Unsigned, NotNull);
	builder->addColumn(ConfigClassInfo::TABLE_NAME, Integer, Unsigned, NotNull);
	builder->addColumn(ConfigStructureVersion::TABLE_NAME, Integer, Unsigned, NotNull);
	builder->addColumn(ConfigVersion::TABLE_NAME, Integer, Unsigned, NotNull);
	builder->addColumn(Site::TABLE_NAME, Integer, Unsigned, NotNull);
	builder->addColumn(Run::TABLE_NAME, Integer, Unsigned, NotNull);
	builder->addColumn(RunCondition::TABLE_NAME, Integer, Unsigned, NotNull);

	createTable(*builder);

	session->once << "INSERT INTO " << dbName << '.' << SEQUENCE_TABLE_NAME <<
		'(' <<
			Config::TABLE_NAME << ',' <<
			ConfigAttribute::TABLE_NAME << ',' <<
			ConfigClass::TABLE_NAME << ',' <<
			ConfigClassInfo::TABLE_NAME << ',' <<
			ConfigStructureVersion::TABLE_NAME << ',' <<
			ConfigVersion::TABLE_NAME << ',' <<
			Site::TABLE_NAME << ',' <<
			Run::TABLE_NAME << ',' <<
			RunCondition::TABLE_NAME <<
		") "
		"VALUES (0,0,0,0,0,0,0,0,0)";
}

/**
 * Creates the table of configurations.
 */
void ConfigDao::createConfigTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(Config::TABLE_NAME);
	builder->addColumn("id", Integer, Unsigned, NotNull);
	builder->addColumn("name", Varchar, SignNA, 255, UNDEFINED, NotNull);
	builder->addPrimaryKey("Pk_Config", "id");
	builder->addUnique("Uk_Config_Name", "name");

	createTable(*builder);
}

/**
 * Creates the table of configuration structure versions.
 */
void ConfigDao::createConfigStructureVersionTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(ConfigStructureVersion::TABLE_NAME);
	builder->addColumn("id", Integer, Unsigned, NotNull);
	builder->addColumn("fkConfig", Integer, Unsigned, NotNull);
	builder->addColumn("creationDate", DateTime, SignNA, NotNull, CurrentDateTime);
	builder->addPrimaryKey("Pk_ConfigStructureVersion", "id");
	builder->addUnique("Uk_ConfigStructureVersion_Date", "fkConfig,creationDate");
	builder->addForeignKey("Fk_StructVer_Config",
			"fkConfig", Config::TABLE_NAME, "id", Cascade, Cascade);

	createTable(*builder);
}

/**
 * Creates the table of configuration versions.
 */
void ConfigDao::createConfigVersionTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(ConfigVersion::TABLE_NAME);
	builder->addColumn("id", Integer, Unsigned, NotNull);
	builder->addColumn("fkStructure", Integer, Unsigned, NotNull);
	builder->addColumn("creationDate", DateTime, SignNA, NotNull, CurrentDateTime);
	builder->addColumn("name", Varchar, SignNA, ConfigVersion::NAME_MAX_LENGTH, UNDEFINED, Null);
	builder->addPrimaryKey("Pk_ConfigVersion", "id");
	builder->addUnique("Uk_ConfigVersion_Date", "fkStructure,creationDate");
	builder->addForeignKey("Fk_ConfigVersion_StructVer",
		"fkStructure", ConfigStructureVersion::TABLE_NAME, "id", Cascade, Cascade);

	createTable(*builder);
}

/**
 * Creates the table of configuration classes.
 */
void ConfigDao::createConfigClassTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(ConfigClass::TABLE_NAME);
	builder->addColumn("id", Integer, Unsigned, NotNull);
	builder->addColumn("fkParent", Integer, Unsigned, Null);
	builder->addColumn("fkCreationStructure", Integer, Unsigned, Null);
	builder->addColumn("name", Varchar, SignNA, ConfigClass::NAME_MAX_LENGTH, UNDEFINED, NotNull);
	builder->addColumn("path", Varchar, SignNA, ConfigClass::PATH_MAX_LENGTH, UNDEFINED, NotNull);
	builder->addPrimaryKey("Pk_ConfigClass", "id");
	builder->addForeignKey("Fk_ConfigClass_Parent",
		"fkParent", ConfigClass::TABLE_NAME, "id", Cascade, Cascade);
	builder->addForeignKey("Fk_ConfigClass_StructVer",
		"fkCreationStructure", ConfigStructureVersion::TABLE_NAME, "id", Cascade, SetNull);

	createTable(*builder);
}

/**
 * Creates the table of configuration class informations.
 */
void ConfigDao::createConfigClassInfoTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(ConfigClassInfo::TABLE_NAME);
	builder->addColumn("id", Integer, Unsigned, NotNull);
	builder->addColumn("fkClass", Integer, Unsigned, NotNull);
	builder->addColumn("fkStructure", Integer, Unsigned, NotNull);
	builder->addColumn("atype", Integer, Unsigned, NotNull);
	builder->addColumn("aunit", Varchar, SignNA, ConfigClassInfo::UNIT_MAX_LENGTH, UNDEFINED, Null);
	builder->addColumn("arange", Varchar, SignNA, ConfigClassInfo::RANGE_MAX_LENGTH, UNDEFINED, Null);
	builder->addPrimaryKey("Pk_ConfigClassInfo", "id");
	builder->addUnique("Uk_ConfigClassInfo", "fkClass,fkStructure");
	builder->addCheck("Ck_ConfigClassInfo_Type", "atype IN (0,1,2,3,4,5)");
	builder->addForeignKey("Fk_ConfigClassInfo_Class",
		"fkClass", ConfigClass::TABLE_NAME, "id", Cascade, Cascade);
	builder->addForeignKey("Fk_ConfigClassInfo_Struct",
			"fkStructure", ConfigStructureVersion::TABLE_NAME, "id", Cascade, Cascade);

	createTable(*builder);
}

/**
 * Creates the table of configuration attributes.
 */
void ConfigDao::createConfigAttributeTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(ConfigAttribute::TABLE_NAME);
	builder->addColumn("id", Integer, Unsigned, NotNull);
	builder->addColumn("fkParent", Integer, Unsigned, Null);
	builder->addColumn("fkVersion", Integer, Unsigned, NotNull);
	builder->addColumn("fkClassInfo", Integer, Unsigned, NotNull);
	builder->addColumn("idx", Varchar, SignNA, ConfigAttribute::INDEX_MAX_LENGTH, UNDEFINED, Null);
	builder->addColumn("value", Varchar, SignNA, ConfigAttribute::VALUE_MAX_LENGTH, UNDEFINED, Null);
	builder->addPrimaryKey("Pk_ConfigAttribute", "id");
	builder->addForeignKey("Fk_ConfigAttribute_Parent",
		"fkParent", ConfigAttribute::TABLE_NAME, "id", Cascade, Cascade);
	builder->addForeignKey("Fk_ConfigAttribute_CfgVer",
		"fkVersion", ConfigVersion::TABLE_NAME, "id", Cascade, Cascade);
	builder->addForeignKey("Fk_ConfigAttribute_ClsInfo",
		"fkClassInfo", ConfigClassInfo::TABLE_NAME, "id", Cascade, Cascade);

	createTable(*builder);
}

/**
 * Creates the table of sites.
 */
void ConfigDao::createSiteTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(Site::TABLE_NAME);
	builder->addColumn("id", Integer, Unsigned, NotNull);
	builder->addColumn("name", Varchar, SignNA, 255, UNDEFINED, NotNull);
	builder->addPrimaryKey("Pk_Site", "id");
	builder->addUnique("Uk_Site_Name", "name");

	createTable(*builder);
}

/**
 * Creates the table of runs.
 */
void ConfigDao::createRunTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(Run::TABLE_NAME);
	builder->addColumn("id", Integer, Unsigned, NotNull);
	builder->addColumn("fkSite", Integer, Unsigned, NotNull);
	builder->addColumn("anumber", Integer, Unsigned, NotNull);
	builder->addPrimaryKey("Pk_Run", "id");
	builder->addUnique("Uk_Run_Number", "fkSite,anumber");
	builder->addForeignKey("Fk_Run_Site", "fkSite", Site::TABLE_NAME, "id", Cascade, Cascade);

	createTable(*builder);
}

/**
 * Creates the table of run conditions.
 */
void ConfigDao::createRunConditionTable()
{
	std::auto_ptr<TableBuilder> builder = createTableBuilder();
	builder->setDatabase(dbName.c_str());
	builder->setTableName(RunCondition::TABLE_NAME);
	builder->addColumn("id", Integer, Unsigned, NotNull);
	builder->addColumn("fkRun", Integer, Unsigned, NotNull);
	builder->addColumn("fkConfigVersion", Integer, Unsigned, NotNull);
	builder->addPrimaryKey("Pk_RunCondition", "id");
	builder->addUnique("Uk_RunConfigVersion", "fkRun,fkConfigVersion");
	builder->addForeignKey("Fk_Run", "fkRun", Run::TABLE_NAME, "id", Cascade, Cascade);
	builder->addForeignKey("Fk_ConfigVersion", "fkConfigVersion", ConfigVersion::TABLE_NAME, "id", Cascade, Cascade);
	createTable(*builder);
}

/**
 * Loads a single ConfigDatabase entity produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns an entity bound with results of the query, or null if it
 * returned an empty result set.
 */
ConfigDatabasePtr ConfigDao::loadDatabase(SOCI::Statement& stmt)
{
	ConfigDatabasePtr database(0);
	std::string name;
	eIndicator indicator = eNoData;

	stmt.exchange(into(name, indicator));
	stmt.defineAndBind();
	stmt.execute();
	stmt.fetch();

	if(indicator == eOK)
	{
		database = new ConfigDatabase;
		database->id = 0;
		database->name = name;
	}
	return database;
}

/**
 * Loads a single Config entity produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns an entity bound with results of the query, or null if it
 * returned an empty result set.
 */
ConfigPtr ConfigDao::loadConfig(SOCI::Statement& stmt)
{
	ConfigPtr config(0);

	UID id;
	std::string name;
	eIndicator indId = eNoData, indName;

	stmt.exchange(into(id, indId));
	stmt.exchange(into(name, indName));
	stmt.defineAndBind();
	stmt.execute();
	stmt.fetch();

	if(indId == eOK)
	{
		config = new Config;
		config->id = id;
		config->name = name;
		config->setDAO(this);
	}
	return config;
}

/**
 * Loads a single ConfigStructureVersion entity produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns an entity bound with results of the query, or null if it
 * returned an empty result set.
 */
ConfigStructureVersionPtr ConfigDao::loadStructure(SOCI::Statement& stmt)
{
	ConfigStructureVersionPtr structVersion(0);
	UID id, fkConfig;
	std::tm creationDate;
	eIndicator indicator = eNoData, dummyInd;

	stmt.exchange(into(id, indicator));
	stmt.exchange(into(fkConfig, dummyInd));
	stmt.exchange(into(creationDate, dummyInd));
	stmt.defineAndBind();
	stmt.execute();
	stmt.fetch();

	if(indicator == eOK)
	{
		structVersion = new ConfigStructureVersion;
		structVersion->id = id;
		structVersion->fkConfig = fkConfig;
		structVersion->creationDate = creationDate;
		structVersion->setDAO(this);
	}
	return structVersion;
}

/**
 * Loads a single ConfigVersion entity produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns an entity bound with results of the query, or null if it
 * returned an empty result set.
 */
ConfigVersionPtr ConfigDao::loadVersion(SOCI::Statement& stmt)
{
	ConfigVersionPtr version(0);

	UID id, fkStructure;
	std::tm creationDate;
	std::string name;
	eIndicator indicator = eNoData, dummyInd;

	stmt.exchange(into(id, indicator));
	stmt.exchange(into(fkStructure, dummyInd));
	stmt.exchange(into(creationDate, dummyInd));
	stmt.exchange(into(name, dummyInd));
	stmt.defineAndBind();
	stmt.execute();
	stmt.fetch();

	if(indicator == eOK)
	{
		version = new ConfigVersion;
		version->id = id;
		version->fkStructure = fkStructure;
		version->creationDate = creationDate;
		version->name = name;
		version->setDAO(this);
	}
	return version;
}

/**
 * Loads a single Site entity produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns an entity bound with results of the query, or null if it
 * returned an empty result set.
 */
SitePtr ConfigDao::loadSite(SOCI::Statement& stmt)
{
	SitePtr site(0);

	UID id;
	std::string name;
	eIndicator indId = eNoData, indName;

	stmt.exchange(into(id, indId));
	stmt.exchange(into(name, indName));
	stmt.defineAndBind();
	stmt.execute();
	stmt.fetch();

	if(indId == eOK)
	{
		site = new Site;
		site->id = id;
		site->name = name;
		site->setDAO(this);
	}
	return site;
}

/**
 * Loads a single Run entity produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns an entity bound with results of the query, or null if it
 * returned an empty result set.
 */
RunPtr ConfigDao::loadRun(SOCI::Statement& stmt)
{
	RunPtr run(0);

	UID id, fkSite, number;
	eIndicator indId = eNoData, dummyInd;

	stmt.exchange(into(id, indId));
	stmt.exchange(into(fkSite, dummyInd));
	stmt.exchange(into(number, dummyInd));
	stmt.defineAndBind();
	stmt.execute();
	stmt.fetch();

	if(indId == eOK)
	{
		run = new Run;
		run->id = id;
		run->fkSite = fkSite;
		run->number = number;
		run->setDAO(this);
	}
	return run;
}

/**
 * Loads a list of Config entities produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns the list of all entities in the resultset.
 */
ConfigPtrSetPtr ConfigDao::loadConfigs(SOCI::Statement& stmt)
{
	ConfigPtrSetPtr configs(new ConfigPtrSet);
	ConfigBatch batch;
	batch.resize(16);

	stmt.exchange(into(batch.getUidBatch(Config::ID)));
	stmt.exchange(into(batch.getStringBatch(Config::NAME)));
	stmt.defineAndBind();
	stmt.execute();

	if(batch.size() != 0)
	{
		while(stmt.fetch())
		{
			for(size_t i = 0; i < batch.size(); ++i)
			{
				ConfigPtr cfg = batch.createEntity(i);
				cfg->setDAO(this);
				configs->insert(cfg);
			}
		}
	}
	return configs;
}

/**
 * Loads a list of ConfigStructureVersion entities produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns the list of all entities in the resultset.
 */
ConfigStructureVersionPtrSetPtr ConfigDao::loadStructures(SOCI::Statement& stmt)
{
	ConfigStructureVersionPtrSetPtr structVersions(new ConfigStructureVersionPtrSet);
	ConfigStructureVersionBatch batch;
	batch.resize(64);

	stmt.exchange(into(batch.getUidBatch(ConfigStructureVersion::ID)));
	stmt.exchange(into(batch.getDateBatch(ConfigStructureVersion::CREATION_DATE)));
	stmt.exchange(into(batch.getUidBatch(ConfigStructureVersion::FK_CONFIG)));
	stmt.defineAndBind();
	stmt.execute();

	if(batch.size() != 0)
	{
		while(stmt.fetch())
		{
			for(size_t i = 0; i < batch.size(); ++i)
			{
				ConfigStructureVersionPtr structVer = batch.createEntity(i);
				structVer->setDAO(this);
				structVersions->insert(structVer);
			}
		}
	}
	return structVersions;
}

/**
 * Loads a list of ConfigVersion entities produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns the list of all entities in the resultset.
 */
ConfigVersionPtrSetPtr ConfigDao::loadVersions(SOCI::Statement& stmt)
{
	ConfigVersionPtrSetPtr versions(new ConfigVersionPtrSet);
	ConfigVersionBatch batch;
	batch.resize(64);
	std::vector<eIndicator> indicName(64);
	StringBatch& nameBatch = batch.getStringBatch(ConfigVersion::NAME);

	stmt.exchange(into(batch.getUidBatch(ConfigVersion::ID)));
	stmt.exchange(into(batch.getDateBatch(ConfigVersion::CREATION_DATE)));
	stmt.exchange(into(nameBatch, indicName));
	stmt.exchange(into(batch.getUidBatch(ConfigVersion::FK_STRUCTURE)));
	stmt.defineAndBind();
	stmt.execute();

	if(batch.size() != 0)
	{
		while(stmt.fetch())
		{
			for(size_t i = 0; i < batch.size(); ++i)
			{
				if(indicName[i] == eNull)
				{
					nameBatch[i].clear();
				}
				ConfigVersionPtr ver = batch.createEntity(i);
				ver->setDAO(this);
				versions->insert(ver);
			}
		}
	}
	return versions;
}

/**
 * Loads a list of ConfigClass entities produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns the list of all entities in the resultset.
 */
ConfigClassPtrSetPtr ConfigDao::loadClasses(SOCI::Statement& stmt)
{
	ConfigClassPtrSetPtr classes(new ConfigClassPtrSet);
	ConfigClassBatch batch;
	batch.resize(64);
	std::vector<eIndicator> indicParent(64);
	std::vector<eIndicator> indicStruct(64);
	UidBatch& fkParentBatch = batch.getUidBatch(ConfigClass::FK_PARENT);
	UidBatch& fkCreatBatch = batch.getUidBatch(ConfigClass::FK_CREATION_STRUCT);

	stmt.exchange(into(batch.getUidBatch(ConfigClass::ID)));
	stmt.exchange(into(fkParentBatch, indicParent));
	stmt.exchange(into(fkCreatBatch, indicStruct));
	stmt.exchange(into(batch.getStringBatch(ConfigClass::NAME)));
	stmt.exchange(into(batch.getStringBatch(ConfigClass::PATH)));
	stmt.defineAndBind();
	stmt.execute();

	if(batch.size() != 0)
	{
		while(stmt.fetch())
		{
			for(size_t i = 0; i < batch.size(); ++i)
			{
				if(indicParent[i] == eNull)
				{
					fkParentBatch[i] = NULL_UID;
				}
				if(indicStruct[i] == eNull)
				{
					fkCreatBatch[i] = NULL_UID;
				}
				ConfigClassPtr cfgClass = batch.createEntity(i);
				cfgClass->setDAO(this);
				classes->insert(cfgClass);
			}
		}
	}
	return classes;
}

/**
 * Loads a list of ConfigClassInfo entities produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns the list of all entities in the resultset.
 */
ConfigClassInfoPtrSetPtr ConfigDao::loadClassInfos(SOCI::Statement& stmt)
{
	ConfigClassInfoPtrSetPtr classInfos(new ConfigClassInfoPtrSet);
	ConfigClassInfoBatch batch;
	batch.resize(64);
	std::vector<eIndicator> indicUnit(64);
	std::vector<eIndicator> indicRange(64);
	StringBatch& unitBatch = batch.getStringBatch(ConfigClassInfo::UNIT);
	StringBatch& rangeBatch = batch.getStringBatch(ConfigClassInfo::RANGE);

	stmt.exchange(into(batch.getUidBatch(ConfigClassInfo::ID)));
	stmt.exchange(into(batch.getUidBatch(ConfigClassInfo::FK_CLASS)));
	stmt.exchange(into(batch.getUidBatch(ConfigClassInfo::FK_STRUCTURE)));
	stmt.exchange(into(batch.getIntBatch(ConfigClassInfo::TYPE)));
	stmt.exchange(into(unitBatch, indicUnit));
	stmt.exchange(into(rangeBatch, indicRange));
	stmt.defineAndBind();
	stmt.execute();

	if(batch.size() != 0)
	{
		while(stmt.fetch())
		{
			for(size_t i = 0; i < batch.size(); ++i)
			{
				if(indicUnit[i] == eNull)
				{
					unitBatch[i].clear();
				}
				if(indicRange[i] == eNull)
				{
					rangeBatch[i].clear();
				}
				ConfigClassInfoPtr info = batch.createEntity(i);
				info->setDAO(this);
				classInfos->insert(info);
			}
		}
	}
	return classInfos;
}

/**
 * Loads a list of ConfigAttribute entities produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns the list of all entities in the resultset.
 */
ConfigAttributePtrSetPtr ConfigDao::loadAttributes(SOCI::Statement& stmt)
{
	ConfigAttributePtrSetPtr attributes(new ConfigAttributePtrSet);
	ConfigAttributeBatch batch;
	batch.resize(256);
	std::vector<eIndicator> indicParent(256);
	std::vector<eIndicator> indicIndex(256);
	std::vector<eIndicator> indicValue(256);

	UidBatch& parentBatch = batch.getUidBatch(ConfigAttribute::FK_PARENT);
	StringBatch& indexBatch = batch.getStringBatch(ConfigAttribute::INDEX);
	StringBatch& valueBatch = batch.getStringBatch(ConfigAttribute::VALUE);

	stmt.exchange(into(batch.getUidBatch(ConfigAttribute::ID)));
	stmt.exchange(into(parentBatch, indicParent));
	stmt.exchange(into(batch.getUidBatch(ConfigAttribute::FK_VERSION)));
	stmt.exchange(into(batch.getUidBatch(ConfigAttribute::FK_CLASSINFO)));
	stmt.exchange(into(indexBatch, indicIndex));
	stmt.exchange(into(valueBatch, indicValue));
	stmt.defineAndBind();
	stmt.execute();

	if(batch.size() != 0)
	{
		while(stmt.fetch())
		{
			for(size_t i = 0; i < batch.size(); ++i)
			{
				if(indicParent[i] == eNull)
				{
					parentBatch[i] = NULL_UID;
				}
				if(indicIndex[i] == eNull)
				{
					indexBatch[i].clear();
				}
				if(indicValue[i] == eNull)
				{
					valueBatch[i].clear();
				}
				ConfigAttributePtr cfgAttr = batch.createEntity(i);
				cfgAttr->setDAO(this);
				attributes->insert(cfgAttr);
			}
		}
	}
	return attributes;
}

/**
 * Loads a list of Site entities produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns the list of all entities in the result set.
 */
SitePtrSetPtr ConfigDao::loadSites(SOCI::Statement& stmt)
{
	SitePtrSetPtr sites(new SitePtrSet);
	SiteBatch batch;
	batch.resize(16);

	stmt.exchange(into(batch.getUidBatch(Site::ID)));
	stmt.exchange(into(batch.getStringBatch(Site::NAME)));
	stmt.defineAndBind();
	stmt.execute();

	if(batch.size() != 0)
	{
		while(stmt.fetch())
		{
			for(size_t i = 0; i < batch.size(); ++i)
			{
				SitePtr site = batch.createEntity(i);
				site->setDAO(this);
				sites->insert(site);
			}
		}
	}
	return sites;
}

/**
 * Loads a list of Run entities produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns the list of all entities in the result set.
 */
RunPtrSetPtr ConfigDao::loadRuns(SOCI::Statement& stmt)
{
	RunPtrSetPtr runs(new RunPtrSet);
	RunBatch batch;
	batch.resize(16);

	stmt.exchange(into(batch.getUidBatch(Run::ID)));
	stmt.exchange(into(batch.getUidBatch(Run::FK_SITE)));
	stmt.exchange(into(batch.getUidBatch(Run::NUMBER)));
	stmt.defineAndBind();
	stmt.execute();

	if(batch.size() != 0)
	{
		while(stmt.fetch())
		{
			for(size_t i = 0; i < batch.size(); ++i)
			{
				RunPtr run = batch.createEntity(i);
				run->setDAO(this);
				runs->insert(run);
			}
		}
	}
	return runs;
}

/**
 * Loads a list of RunCondition entities produced by a prepared statement.
 * @param stmt a statement initialized with a query and optionally with input parameters.
 * @return Returns the list of all entities in the result set.
 */
RunConditionPtrSetPtr ConfigDao::loadRunConditions(SOCI::Statement& stmt)
{
	RunConditionPtrSetPtr conditions(new RunConditionPtrSet);
	RunConditionBatch batch;
	batch.resize(16);

	stmt.exchange(into(batch.getUidBatch(RunCondition::ID)));
	stmt.exchange(into(batch.getUidBatch(RunCondition::FK_RUN)));
	stmt.exchange(into(batch.getUidBatch(RunCondition::FK_CONFIGVERSION)));
	stmt.defineAndBind();
	stmt.execute();

	if(batch.size() != 0)
	{
		while(stmt.fetch())
		{
			for(size_t i = 0; i < batch.size(); ++i)
			{
				RunConditionPtr condition = batch.createEntity(i);
				condition->setDAO(this);
				conditions->insert(condition);
			}
		}
	}
	return conditions;
}

/**
 * Returns the list of configurations in a database.
 */
ConfigPtrSetPtr ConfigDao::getConfigs()
{
	Statement stmt = (session->prepare <<
		"SELECT id,name "
		"FROM " << dbName << '.' << Config::TABLE_NAME);
	return loadConfigs(stmt);
}

/**
 * Returns the list of structure versions of a configuration.
 */
ConfigStructureVersionPtrSetPtr ConfigDao::getStructureVersions(Config& cfg)
{
	Statement stmt = (session->prepare <<
		"SELECT id,creationDate,fkConfig "
		"FROM " << dbName << '.' << ConfigStructureVersion::TABLE_NAME << " "
		"WHERE fkConfig=" << cfg.id);
	return loadStructures(stmt);
}

/**
 * Returns the list of versions of a configuration.
 */
ConfigVersionPtrSetPtr ConfigDao::getVersionsOfConfig(Config& cfg)
{
	Statement stmt = (session->prepare <<
		"SELECT v.id,v.creationDate,v.name,v.fkStructure "
		"FROM " << dbName << '.' << ConfigVersion::TABLE_NAME << " v "
			"INNER JOIN " <<
			dbName << '.' << ConfigStructureVersion::TABLE_NAME << " s "
			"ON v.fkStructure=s.id "
		"WHERE s.fkConfig=" << cfg.id);
	return loadVersions(stmt);
}

/**
 * Returns the list of configuration version sharing the same structure.
 */
ConfigVersionPtrSetPtr ConfigDao::getVersionsOfStructure(ConfigStructureVersion& structure)
{
	Statement stmt = (session->prepare <<
		"SELECT id,creationDate,name,fkStructure "
		"FROM " << dbName << '.' << ConfigVersion::TABLE_NAME << " "
		"WHERE fkStructure=" << structure.id);
	return loadVersions(stmt);
}

/**
 * Returns the list of configuration classes in use for a given structure version.
 */
ConfigClassPtrSetPtr ConfigDao::getClassesOfStructure(ConfigStructureVersion& structure)
{
	Statement stmt = (session->prepare <<
		"SELECT c.id,c.fkParent,c.fkCreationStructure,c.name,c.path "
		"FROM " << dbName << '.' << ConfigClassInfo::TABLE_NAME << " i "
			"INNER JOIN " << dbName << '.' << ConfigClass::TABLE_NAME << " c "
			"ON i.fkClass=c.id "
		"WHERE i.fkStructure=" << structure.id);
	return loadClasses(stmt);
}

/**
 * Returns the list of configuration classes related to a configuration.
 */
ConfigClassPtrSetPtr ConfigDao::getClassesOfConfig(Config& cfg)
{
	Statement stmt = (session->prepare <<
		"SELECT c.id,c.fkParent,c.fkCreationStructure,c.name,c.path "
		"FROM " << dbName << '.' << ConfigClass::TABLE_NAME << " c "
			"INNER JOIN " << dbName << '.' << ConfigStructureVersion::TABLE_NAME << " s "
			"ON c.fkCreationStructure=s.id "
		"WHERE s.fkConfig=" << cfg.id);
	return loadClasses(stmt);
}

/**
 * Returns the list of configuration class informations for a given structure version.
 */
ConfigClassInfoPtrSetPtr ConfigDao::getClassInfosOfStructure(ConfigStructureVersion& structure)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkClass,fkStructure,atype,aunit,arange "
		"FROM " << dbName << '.' << ConfigClassInfo::TABLE_NAME << " "
		"WHERE fkStructure=" << structure.id);
	return loadClassInfos(stmt);
}

/**
 * Returns the list of the configuration attributes belonging to a version.
 */
ConfigAttributePtrSetPtr ConfigDao::getAttributesOfVersion(ConfigVersion& version)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkParent,fkVersion,fkClassInfo,idx,value "
		"FROM " << dbName << '.' << ConfigAttribute::TABLE_NAME << " "
		"WHERE fkVersion=" << version.id);
	return loadAttributes(stmt);
}

/**
 * Returns the list of the runs belonging to a site.
 */
RunPtrSetPtr ConfigDao::getRunsOfSite(Site& site)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkSite,anumber "
		"FROM " << dbName << '.' << Run::TABLE_NAME << " "
		"WHERE fkSite=" << site.id);
	return loadRuns(stmt);
}

/**
 * Returns the list of run conditions belonging to a specified run.
 */
RunConditionPtrSetPtr ConfigDao::getConditionsOfRun(Run& run)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkRun,fkConfigVersion "
		"FROM " << dbName << '.' << RunCondition::TABLE_NAME << " "
		"WHERE fkRun=" << run.id);
	return loadRunConditions(stmt);
}

/**
 * Returns the list of versions of a configuration related to a run.
 */
ConfigVersionPtrSetPtr ConfigDao::getVersionsOfRun(Run& run)
{
	Statement stmt = (session->prepare <<
		"SELECT v.id,v.creationDate,v.name,v.fkStructure "
		"FROM " << dbName << '.' << ConfigVersion::TABLE_NAME << " v "
			"INNER JOIN " <<
			dbName << '.' << RunCondition::TABLE_NAME << " c "
			"ON v.id=c.fkConfigVersion "
		"WHERE c.fkRun=" << run.id);
	return loadVersions(stmt);
}

/**
 * Returns the list of sites in a database.
 */
SitePtrSetPtr ConfigDao::getSites()
{
	Statement stmt = (session->prepare <<
		"SELECT id,name "
		"FROM " << dbName << '.' << Site::TABLE_NAME);
	return loadSites(stmt);
}

/**
 * Loads the Config entity with the specified id from the database.
 * @return Returns the entity or null if an entity with such an id does not exist.
 */
ConfigPtr ConfigDao::getConfig(UID id)
{
	Statement stmt = (session->prepare <<
		"SELECT id,name "
		"FROM " << dbName << '.' << Config::TABLE_NAME << " "
		"WHERE id=" << id);
	return loadConfig(stmt);
}

/**
 * Loads the ConfigStructureVersion entity with the specified id from the database.
 * @return Returns the entity or null if an entity with such an id does not exist.
 */
ConfigStructureVersionPtr ConfigDao::getStructure(UID id)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkConfig,creationDate "
		"FROM " << dbName << '.' << ConfigStructureVersion::TABLE_NAME << " "
		"WHERE id=" << id);
	return loadStructure(stmt);
}

/**
 * Loads the ConfigVersion entity with the specified id from the database.
 * @return Returns the entity or null if an entity with such an id does not exist.
 */
ConfigVersionPtr ConfigDao::getVersion(UID id)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkStructure,creationDate,name "
		"FROM " << dbName << '.' << ConfigVersion::TABLE_NAME << " "
		"WHERE id=" << id);
	return loadVersion(stmt);
}

/**
 * Returns the configuration bearing the specified name, if it exists.
 */
ConfigPtr ConfigDao::searchConfig(const std::string& cfgName)
{
	Statement stmt = (session->prepare <<
		"SELECT id,name "
		"FROM " << dbName << '.' << Config::TABLE_NAME << " "
		"WHERE name=:cfgName",
		use(const_cast<std::string&>(cfgName)));
	return loadConfig(stmt);
}

/**
 * Returns the most recent structure version of a configuration.
 */
ConfigStructureVersionPtr ConfigDao::getLastStructureVersion(Config& cfg)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkConfig,creationDate "
		"FROM " << dbName << '.' << ConfigStructureVersion::TABLE_NAME << " "
		"WHERE creationDate=("
			"SELECT MAX(creationDate) "
			"FROM " << dbName << '.' << ConfigStructureVersion::TABLE_NAME << " "
			"WHERE fkConfig=" << cfg.id << ')');
	return loadStructure(stmt);
}

/**
 * Returns the structure version matching exactly the specified date.
 */
ConfigStructureVersionPtr ConfigDao::getStructureVersionAt(const std::tm& date)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkConfig,creationDate "
		"FROM " << dbName << '.' << ConfigStructureVersion::TABLE_NAME << " "
		"WHERE creationDate=:creationDate",
		use(const_cast<std::tm&>(date)));
	return loadStructure(stmt);
}

/**
 * Returns the most recent version of a configuration.
 */
ConfigVersionPtr ConfigDao::getLastVersion(Config& cfg)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkStructure,creationDate,name "
		"FROM " << dbName << '.' << ConfigVersion::TABLE_NAME << " "
		"WHERE creationDate=("
			"SELECT MAX(v.creationDate) "
			"FROM " << dbName << '.' << ConfigVersion::TABLE_NAME << " v "
				"INNER JOIN " << dbName << '.' << ConfigStructureVersion::TABLE_NAME << " s "
				"ON v.fkStructure=s.id "
			"WHERE s.fkConfig=" << cfg.id << ')');
	return loadVersion(stmt);
}

/**
 * Returns the version record matching exactly the specified date.
 */
ConfigVersionPtr ConfigDao::getVersionAt(const std::tm& date)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkStructure,creationDate,name "
		"FROM " << dbName << '.' << ConfigVersion::TABLE_NAME << " "
		"WHERE creationDate=:creationDate",
		use(const_cast<std::tm&>(date)));
	return loadVersion(stmt);
}

/**
 * Loads the Site entity with the specified id from the database.
 * @return Returns the entity or null if an entity with such an id does not exist.
 */
SitePtr ConfigDao::getSite(UID id)
{
	Statement stmt = (session->prepare <<
		"SELECT id,name "
		"FROM " << dbName << '.' << Site::TABLE_NAME << " "
		"WHERE id=" << id);
	return loadSite(stmt);
}

/**
 * Loads the site bearing the specified name, if it exists.
 * @return Returns the entity or null if an entity with such a name does not exist.
 */
SitePtr ConfigDao::searchSite(const std::string& siteName)
{
	Statement stmt = (session->prepare <<
		"SELECT id,name "
		"FROM " << dbName << '.' << Site::TABLE_NAME << " "
		"WHERE name=:siteName",
		use(const_cast<std::string&>(siteName)));
	return loadSite(stmt);
}

/**
 * Loads the Run entity with the specified id from the database.
 * @return Returns the entity or null if an entity with such an id does not exist.
 */
RunPtr ConfigDao::getRun(UID id)
{
	Statement stmt = (session->prepare <<
		"SELECT id,fkSite,anumber "
		"FROM " << dbName << '.' << Run::TABLE_NAME << " "
		"WHERE id=" << id);
	return loadRun(stmt);
}

/**
 * Loads the Run entity with the specified site name and run number from the database.
 * @return Returns the entity or null if an entity with such a (site name, run number) does not exist.
 */
RunPtr ConfigDao::searchRun(const std::string& siteName, const UID& runNumber)
{
	Statement stmt = (session->prepare <<
		"SELECT r.id,r.fkSite,r.anumber "
		"FROM " << dbName << '.' << Run::TABLE_NAME << " r "
			"INNER JOIN " << dbName << '.' << Site::TABLE_NAME << " s "
			"ON r.fkSite=s.id "
		"WHERE s.name=:siteName AND r.anumber=" << runNumber,
		use(const_cast<std::string&>(siteName)));
	return loadRun(stmt);
}

/**
 * Allocates primary key numbers to identity new entities that will be inserted.
 * @param tableName the table where reservations shall be made.
 * @param count the number of entities.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRange(const char* tableName, size_t count)
{
	UID rangeBegin;
	Transaction transac(*session);
	session->once <<
		"SELECT " << tableName << " "
		"FROM " << dbName << '.' << SEQUENCE_TABLE_NAME << " FOR UPDATE",
		into(rangeBegin);
	session->once <<
		"UPDATE " << dbName << '.' << SEQUENCE_TABLE_NAME << " "
		"SET " << tableName << '=' << tableName << '+' << count;
	transac.commit();
	return rangeBegin;
}

/**
 * Allocates primary key numbers to identity new configurations that will be inserted.
 * @param count the number of configurations.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRangeForConfigs(size_t count)
{
	return allocIdRange(Config::TABLE_NAME, count);
}

/**
 * Allocates primary key numbers to identity new structure versions that will be inserted.
 * @param count the number of structure versions.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRangeForStructureVersions(size_t count)
{
	return allocIdRange(ConfigStructureVersion::TABLE_NAME, count);
}

/**
 * Allocates primary key numbers to identity new versions that will be inserted.
 * @param count the number of structure versions.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRangeForVersions(size_t count)
{
	return allocIdRange(ConfigVersion::TABLE_NAME, count);
}

/**
 * Allocates primary key numbers to identity new classes that will be inserted.
 * @param count the number of classes.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRangeForClasses(size_t count)
{
	return allocIdRange(ConfigClass::TABLE_NAME, count);
}

/**
 * Allocates primary key numbers to identity new class infos that will be inserted.
 * @param count the number of class infos.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRangeForClassInfos(size_t count)
{
	return allocIdRange(ConfigClassInfo::TABLE_NAME, count);
}

/**
 * Allocates primary key numbers to identity new attributes that will be inserted.
 * @param count the number of attributes.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRangeForAttributes(size_t count)
{
	return allocIdRange(ConfigAttribute::TABLE_NAME, count);
}

/**
 * Allocates primary key numbers to identity new sites that will be inserted.
 * @param count the number of sites.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRangeForSites(size_t count)
{
	return allocIdRange(Site::TABLE_NAME, count);
}

/**
 * Allocates primary key numbers to identity new runs that will be inserted.
 * @param count the number of runs.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRangeForRuns(size_t count)
{
	return allocIdRange(Run::TABLE_NAME, count);
}

/**
 * Allocates primary key numbers to identity new run conditions that will be inserted.
 * @param count the number of run conditions.
 * @return Returns the lower bound of the usable range ]id, id + count] (so it must
 * be pre-incremented to get the first usable id).
 */
UID ConfigDao::allocIdRangeForRunConditions(size_t count)
{
	return allocIdRange(RunCondition::TABLE_NAME, count);
}

/**
 * Allocates primary key numbers for all entity types in a single call.
 * @param counts in an array containing the number of ids to alloc for each table.
 * Use the enum Tables to index this array. Zero is a perfectly valid number if
 * no id should be allocated for some of the tables. In this case, the returned
 * id for this table should not be used to initialize an entity.
 * @param ids out an array that will be filled with the first id of the range for
 * each table.
 */
void ConfigDao::allocIdRanges(const UID (&counts)[TableCount], UID (&ids)[TableCount])
{
	Transaction transac(*session);
	session->once <<
		"SELECT " <<
			TABLES[TConfig] << ',' <<
			TABLES[TStructureVersion] << ',' <<
			TABLES[TVersion] << ',' <<
			TABLES[TClassInfo] << ',' <<
			TABLES[TClass] << ',' <<
			TABLES[TAttribute] << ',' <<
			TABLES[TSite] << ',' <<
			TABLES[TRun] << ',' <<
			TABLES[TRunCondition] << " "
		"FROM " << dbName << '.' << SEQUENCE_TABLE_NAME << " FOR UPDATE",
		into(ids[TConfig]),
		into(ids[TStructureVersion]),
		into(ids[TVersion]),
		into(ids[TClassInfo]),
		into(ids[TClass]),
		into(ids[TAttribute]),
		into(ids[TSite]),
		into(ids[TRun]),
		into(ids[TRunCondition]);
	session->once <<
		"UPDATE " << dbName << '.' << SEQUENCE_TABLE_NAME << " "
		"SET " <<
			TABLES[TConfig] << '=' << TABLES[TConfig] << '+' << counts[TConfig] << ',' <<
			TABLES[TStructureVersion] << '=' << TABLES[TStructureVersion] << '+' << counts[TStructureVersion] << ',' <<
			TABLES[TVersion] << '=' << TABLES[TVersion] << '+' << counts[TVersion] << ',' <<
			TABLES[TClassInfo] << '=' << TABLES[TClassInfo] << '+' << counts[TClassInfo] << ',' <<
			TABLES[TClass] << '=' << TABLES[TClass] << '+' << counts[TClass] << ',' <<
			TABLES[TAttribute] << '=' << TABLES[TAttribute] << '+' << counts[TAttribute] << ',' <<
			TABLES[TSite] << '=' << TABLES[TSite] << '+' << counts[TSite] << ',' <<
			TABLES[TRun] << '=' << TABLES[TRun] << '+' << counts[TRun] << ',' <<
			TABLES[TRunCondition] << '=' << TABLES[TRunCondition] << '+' << counts[TRunCondition];
	transac.commit();
}

/**
 * Inserts a configuration into the database.
 */
void ConfigDao::insert(Config& config)
{
	session->once <<
		"INSERT INTO " << dbName << '.' << Config::TABLE_NAME << "(id,name) "
		"VALUES(NULLIF(:id,0),:name)",
		use(config.id),
		use(config.name);
}

/**
 * Inserts a structure version into the database.
 */
void ConfigDao::insert(ConfigStructureVersion& structVersion)
{
	session->once <<
		"INSERT INTO " << dbName << '.' << ConfigStructureVersion::TABLE_NAME <<
			"(id,creationDate,fkConfig) "
		"VALUES(NULLIF(:id,0),:creationDate,NULLIF(:fkConfig,0))",
		use(structVersion.id),
		use(structVersion.creationDate),
		use(structVersion.fkConfig);
}

/**
 * Inserts a version into the database.
 */
void ConfigDao::insert(ConfigVersion& version)
{
	session->once <<
	    "INSERT INTO " << dbName << '.' << ConfigVersion::TABLE_NAME <<
	    	"(id,creationDate,name,fkStructure) "
	    "VALUES(NULLIF(:id,0),:creationDate,:name,NULLIF(:fkStructure,0))",
        use(version.id),
        use(version.creationDate),
        use(version.name),
        use(version.fkStructure);
}

/**
 * Inserts a site into the database.
 */
void ConfigDao::insert(Site& site)
{
	session->once <<
		"INSERT INTO " << dbName << '.' << Site::TABLE_NAME << "(id,name) "
		"VALUES(NULLIF(:id,0),:name)",
		use(site.id),
		use(site.name);
}

/**
 * Inserts a run into the database.
 */
void ConfigDao::insert(Run& run)
{
	Transaction transac(*session);
	// Get maximum run number for site
	UID maxRunNumber=-1;
	eIndicator ind(eNoData);
	session->once <<
		"SELECT MAX(anumber) FROM " << dbName << '.' << Run::TABLE_NAME << " WHERE fkSite=:siteId",
		into(maxRunNumber, ind),
		use(run.fkSite);

	// Set run number to next available number
	run.number = 0;
	if (ind == eOK)
	{
			run.number = maxRunNumber+1;
	}

	// Insert run into database
	session->once <<
		"INSERT INTO " << dbName << '.' << Run::TABLE_NAME << "(id,fkSite,anumber) "
		"VALUES(NULLIF(:id,0),:fkSite,:anumber)",
		use(run.id),
		use(run.fkSite),
		use(run.number);
	transac.commit();
}

/**
 * Inserts a run condition (e.g. an association between a run and a version of a configuration) into the database.
 */
void ConfigDao::insert(RunCondition& condition)
{
	session->once <<
	    "INSERT INTO " << dbName << '.' << RunCondition::TABLE_NAME <<
	    	"(id,fkRun,fkConfigVersion) "
	    "VALUES(NULLIF(:id,0),NULLIF(:fkRun,0),NULLIF(:fkConfigVersion,0))",
        use(condition.id),
        use(condition.fkRun),
        use(condition.fkConfigVersion);
}

/**
 * Inserts a batch of configuration into the database.
 */
void ConfigDao::insert(ConfigBatch& batch)
{
	if(batch.size() > 0)
	{
		Statement st = (session->prepare <<
	        "INSERT INTO " << dbName << '.' << Config::TABLE_NAME << "(id,name) "
	        "VALUES(NULLIF(:id,0),:name)",
	        use(batch.getUidBatch(Config::ID)),
	        use(batch.getStringBatch(Config::NAME)));
		st.execute(true);
	}
}

/**
 * Inserts a batch of structure versions into the database.
 */
void ConfigDao::insert(ConfigStructureVersionBatch& batch)
{
	if(batch.size() > 0)
	{
		Statement st = (session->prepare <<
	        "INSERT INTO " << dbName << '.' << ConfigStructureVersion::TABLE_NAME <<
	        	"(id,creationDate,fkConfig) "
	        "VALUES(NULLIF(:id,0),:creationDate,NULLIF(:fkConfig,0))",
	        use(batch.getUidBatch(ConfigStructureVersion::ID)),
	        use(batch.getDateBatch(ConfigStructureVersion::CREATION_DATE)),
	        use(batch.getUidBatch(ConfigStructureVersion::FK_CONFIG)));
		st.execute(true);
	}
}

/**
 * Inserts a batch of versions into the database.
 */
void ConfigDao::insert(ConfigVersionBatch& batch)
{
	if(batch.size() > 0)
	{
		Statement st = (session->prepare <<
	        "INSERT INTO " << dbName << '.' << ConfigVersion::TABLE_NAME <<
	        	"(id,creationDate,name,fkStructure) "
	        "VALUES(NULLIF(:id,0),:creationDate,:name,NULLIF(:fkStructure,0))",
	        use(batch.getUidBatch(ConfigVersion::ID)),
	        use(batch.getDateBatch(ConfigVersion::CREATION_DATE)),
	        use(batch.getStringBatch(ConfigVersion::NAME)),
	        use(batch.getUidBatch(ConfigVersion::FK_STRUCTURE)));
		st.execute(true);
	}
}

/**
 * Inserts a batch of classes into the database.
 */
void ConfigDao::insert(ConfigClassBatch& batch)
{
	if(batch.size() > 0)
	{
		Statement st = (session->prepare <<
	        "INSERT INTO " << dbName << '.' << ConfigClass::TABLE_NAME <<
	        	"(id,fkParent,fkCreationStructure,name,path) "
	        "VALUES(NULLIF(:id,0),NULLIF(:fkParent,0),NULLIF(:fkCreationStructure,0),:name,:path)",
	        use(batch.getUidBatch(ConfigClass::ID)),
	        use(batch.getUidBatch(ConfigClass::FK_PARENT)),
	        use(batch.getUidBatch(ConfigClass::FK_CREATION_STRUCT)),
	        use(batch.getStringBatch(ConfigClass::NAME)),
	        use(batch.getStringBatch(ConfigClass::PATH)));
		st.execute(true);
	}
}

/**
 * Inserts a batch of class infos into the database.
 */
void ConfigDao::insert(ConfigClassInfoBatch& batch)
{
	if(batch.size() > 0)
	{
		Statement st = (session->prepare <<
	        "INSERT INTO " << dbName << '.' << ConfigClassInfo::TABLE_NAME <<
	        	"(id,fkClass,fkStructure,atype,aunit,arange) "
	        "VALUES(NULLIF(:id,0),NULLIF(:fkClass,0),NULLIF(:fkStructure,0),:type,:unit,:range)",
	        use(batch.getUidBatch(ConfigClassInfo::ID)),
	        use(batch.getUidBatch(ConfigClassInfo::FK_CLASS)),
	        use(batch.getUidBatch(ConfigClassInfo::FK_STRUCTURE)),
	        use(batch.getIntBatch(ConfigClassInfo::TYPE)),
	        use(batch.getStringBatch(ConfigClassInfo::UNIT)),
	        use(batch.getStringBatch(ConfigClassInfo::RANGE)));
		st.execute(true);
	}
}

/**
 * Inserts a batch of attributes into the database.
 */
void ConfigDao::insert(ConfigAttributeBatch& batch)
{
	if(batch.size() > 0)
	{
		Statement st = (session->prepare <<
	        "INSERT INTO " << dbName << '.' << ConfigAttribute::TABLE_NAME <<
	        	"(id,fkParent,fkVersion,fkClassInfo,idx,value) "
	        "VALUES(NULLIF(:id,0),NULLIF(:fkParent,0),NULLIF(:fkVersion,0),"
	        	"NULLIF(:fkClassInfo,0),:idx,:value)",
	        use(batch.getUidBatch(ConfigAttribute::ID)),
	        use(batch.getUidBatch(ConfigAttribute::FK_PARENT)),
	        use(batch.getUidBatch(ConfigAttribute::FK_VERSION)),
	        use(batch.getUidBatch(ConfigAttribute::FK_CLASSINFO)),
	        use(batch.getStringBatch(ConfigAttribute::INDEX)),
	        use(batch.getStringBatch(ConfigAttribute::VALUE)));
		st.execute(true);
	}
}

/**
 * Inserts a batch of sites into the database.
 */
void ConfigDao::insert(SiteBatch& batch)
{
	if(batch.size() > 0)
	{
		Statement st = (session->prepare <<
	        "INSERT INTO " << dbName << '.' << Site::TABLE_NAME <<
	        	"(id,name) "
	        "VALUES(NULLIF(:id,0),:name)",
	        use(batch.getUidBatch(Site::ID)),
	        use(batch.getStringBatch(Site::NAME)));
		st.execute(true);
	}
}

/**
 * Inserts an entire configuration loaded in entity batches, using a transaction.
 */
void ConfigDao::insertFullConfig(ConfigStructureVersionBatch& sversions,
			ConfigVersionBatch& versions,
			ConfigClassBatch& classes,
			ConfigClassInfoBatch& classInfos,
			ConfigAttributeBatch& attributes)
{
	Transaction transac(*session);
	insert(sversions);
	insert(versions);
	insert(classes);
	insert(classInfos);
	insert(attributes);
	transac.commit();
}

/**
 * Removes a configuration and all its associated records in cascade.
 */
void ConfigDao::remove(Config& config)
{
	session->once <<
		"DELETE FROM " << dbName << '.' << Config::TABLE_NAME << " "
		"WHERE id=" << config.id;
}

/**
 * Removes a structure version and all its associated records in cascade.
 */
void ConfigDao::remove(ConfigStructureVersion& structVersion)
{
	session->once <<
		"DELETE FROM " << dbName << '.' << ConfigStructureVersion::TABLE_NAME << " "
		"WHERE id=" << structVersion.id;
}

/**
 * Removes a version and all its config attributes in cascade.
 * If the parent structure is not referenced any more, it is also removed.
 */
void ConfigDao::remove(ConfigVersion& version)
{
	Transaction transac(*session);
	session->once <<
		"DELETE FROM " << dbName << '.' << ConfigVersion::TABLE_NAME << " "
		"WHERE id=" << version.id;
	transac.commit();
}

/**
 * Removes a site and all its runs and their run conditions in cascade.
 */
void ConfigDao::remove(Site& site)
{
	session->once <<
		"DELETE FROM " << dbName << '.' << Site::TABLE_NAME << " "
		"WHERE id=" << site.id;
}

/**
 * Removes a run and all its run conditions in cascade.
 */
void ConfigDao::remove(Run& run)
{
	session->once <<
		"DELETE FROM " << dbName << '.' << Run::TABLE_NAME << " "
		"WHERE id=" << run.id;
}

/**
 * Removes a run condition.
 */
void ConfigDao::remove(RunCondition& condition)
{
	session->once <<
		"DELETE FROM " << dbName << '.' << RunCondition::TABLE_NAME << " "
		"WHERE id=" << condition.id;
}

/**
 * Changes the name of the specified configuration record.
 */
void ConfigDao::rename(Config& config, const std::string& newName)
{
	session->once <<
		"UPDATE " << dbName << '.' << Config::TABLE_NAME << " "
		"SET name=:name "
		"WHERE id=:id",
		use(const_cast<std::string&>(newName)),
		use(config.id);
	config.name = newName;
}

/**
 * Changes the name of the specified version record.
 */
void ConfigDao::rename(ConfigVersion& version, const std::string& newName)
{
	session->once <<
		"UPDATE " << dbName << '.' << ConfigVersion::TABLE_NAME << " "
		"SET name=:name "
		"WHERE id=:id",
		use(const_cast<std::string&>(newName)),
		use(version.id);
	version.name = newName;
}

/**
 * Changes the name of the specified site record.
 */
void ConfigDao::rename(Site& site, const std::string& newName)
{
	session->once <<
		"UPDATE " << dbName << '.' << Site::TABLE_NAME << " "
		"SET name=:name "
		"WHERE id=:id",
		use(const_cast<std::string&>(newName)),
		use(site.id);
	site.name = newName;
}

/**
 * @fn std::auto_ptr<TableBuilder> ConfigDao::createTableBuilder()
 * Creates and return an object to build CREATE TABLE queries in a RDMS neutral way.
 */

/**
 * Executes queries of the table builder to create its table and all related
 * utility objects.
 */
void ConfigDao::createTable(TableBuilder& builder)
{
	session->once << builder.buildQuery();

	const StringBatch& sideQueries = builder.getSideQueries();
	for(size_t i = 0; i < sideQueries.size(); ++i)
	{
		session->once << sideQueries[i];
	}
}

/**
 * Prints in a stream the list of privilege from a combination of flags.
 */
std::string ConfigDao::formatPrivileges(int privilege)
{
	std::string privs;
	int itemCount = 0;
	if(privilege & Select)
	{
		privs += "SELECT";
		++itemCount;
	}
	if(privilege & Insert)
	{
		if(itemCount++ > 0) privs += ',';
		privs += "INSERT";
	}
	if(privilege & Update)
	{
		if(itemCount++ > 0) privs += ',';
		privs += "UPDATE";
	}
	if(privilege & Delete)
	{
		if(itemCount++ > 0)	privs += ',';
		privs += "DELETE";
	}
	return privs;
}

}
}
}
