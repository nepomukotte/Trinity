#ifndef CCFG_IO_DB_CONFIGDAO_H_INCLUDED
#define CCFG_IO_DB_CONFIGDAO_H_INCLUDED

#include "Config.h"
#include "ConfigAttribute.h"
#include "ConfigClass.h"
#include "ConfigClassInfo.h"
#include "ConfigDatabase.h"
#include "ConfigStructureVersion.h"
#include "ConfigVersion.h"
#include "Site.h"
#include "Run.h"
#include "RunCondition.h"
#include <string>
#include <memory>

// forward declarations
namespace SOCI { class Session; class Statement; }

namespace CCfg
{
namespace Io
{
namespace Db
{
// enums
enum ConfigUserRole {Expert = 0, Shifter, Analyst, NoAccess};
enum PrivilegeFlags {Select = 0x1, Insert = 0x2, Update = 0x4, Delete = 0x8};
enum TableId {TConfig = 0, TStructureVersion, TVersion, TClassInfo, TClass, TAttribute, TSite, TRun, TRunCondition, TableCount};
#define ROLE_COUNT	Analyst+1

// forward declarations
class TableBuilder;
class ConfigBatch;
class ConfigAttributeBatch;
class ConfigClassBatch;
class ConfigClassInfoBatch;
class ConfigStructureVersionBatch;
class ConfigVersionBatch;
class SiteBatch;


////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is a Data Access Object: it encapsulates every database query in an interface
 * implemented for each RDMS.
 *
 * @author fchateau
 */
class ConfigDao
{
protected:
	std::auto_ptr<SOCI::Session> session;
	std::string dbName;

public:

	static const char SEQUENCE_TABLE_NAME[];
	static const char* TABLES[TableCount+1];
	static const unsigned TABLE_PRIVILEGES[TableCount+1][ROLE_COUNT];

	ConfigDao();
	virtual ~ConfigDao();

	// @{ Connection management
	virtual void connect(const std::string& connectionString) = 0;
	bool isConnected() const { return session.get() != 0; }
	void setWorkingDatabase(const std::string& dbName) { this->dbName = dbName; }
	const std::string& getWorkingDatabase() const { return dbName; }
	// @}

	// @{ Database and user management
	virtual ConfigDatabasePtr getCurrentDatabase() = 0;
	virtual ConfigDatabasePtr searchDatabase(const std::string& dbName) = 0;
	virtual ConfigDatabasePtrSetPtr getDatabases() = 0;
	virtual std::vector<std::string> getUserNames() = 0;
	virtual void createDatabaseAndUser(const std::string& name, const std::string& pwd) = 0;
	virtual void removeDatabaseAndUser(const std::string& name) = 0;
	virtual void changeUserPassword(const std::string& userName, const std::string& newPassword) = 0;
	// @}

	// @{ Privileges management
	virtual int getPrivileges(const std::string& userName, const std::string& tableName) = 0;
	void grantPrivileges(const std::string& userName, const std::string& tableName, int privileges);
	void revokePrivileges(const std::string& userName, int privileges);
	void revokePrivileges(const std::string& userName, const std::string& tableName, int privileges);
	virtual void revokePrivileges(const std::string& userName) = 0;
	// @}

	// @{ Roles management
	ConfigUserRole getRole(const std::string& userName);
	ConfigUserRole getRole(const std::string& userName, int tableIndex);
	void setRole(const std::string& userName, ConfigUserRole role);
	void grantRole(const std::string& userName, ConfigUserRole role);
	virtual bool isSchemaAdmin(const std::string& userName) = 0;
	// @}

	// @{ Tables management
	virtual bool checkTableExist(const char* tableName) = 0;
	bool checkHasAllTables();
	void dropAllTables(bool nothrow = false);
	void dropTable(const char* tableName, bool nothrow = false);
	void replaceAllTables();
	void createAllTables();
	void createConfigSequenceTable();
	void createConfigTable();
	void createConfigStructureVersionTable();
	void createConfigVersionTable();
	void createConfigClassTable();
	void createConfigClassInfoTable();
	void createConfigAttributeTable();
	void createSiteTable();
	void createRunTable();
	void createRunConditionTable();
	// @}

	// @{ entity set retrieval
	ConfigPtrSetPtr getConfigs();
	ConfigStructureVersionPtrSetPtr getStructureVersions(Config& cfg);
	ConfigVersionPtrSetPtr getVersionsOfConfig(Config& cfg);
	ConfigVersionPtrSetPtr getVersionsOfStructure(ConfigStructureVersion& structure);
	ConfigClassPtrSetPtr getClassesOfStructure(ConfigStructureVersion& structure);
	ConfigClassPtrSetPtr getClassesOfConfig(Config& cfg);
	ConfigClassInfoPtrSetPtr getClassInfosOfStructure(ConfigStructureVersion& structure);
	ConfigAttributePtrSetPtr getAttributesOfVersion(ConfigVersion& version);
	SitePtrSetPtr getSites();
	RunPtrSetPtr getRunsOfSite(Site& site);
	RunConditionPtrSetPtr getConditionsOfRun(Run& run);
	ConfigVersionPtrSetPtr getVersionsOfRun(Run& run);
	// @}

	// @{ individual entity retrieval
	ConfigPtr getConfig(UID id);
	ConfigStructureVersionPtr getStructure(UID id);
	ConfigVersionPtr getVersion(UID id);
	ConfigPtr searchConfig(const std::string& cfgName);
	ConfigStructureVersionPtr getLastStructureVersion(Config& cfg);
	ConfigStructureVersionPtr getStructureVersionAt(const std::tm& date);
	ConfigVersionPtr getLastVersion(Config& cfg);
	ConfigVersionPtr getVersionAt(const std::tm& date);
	SitePtr getSite(UID id);
	SitePtr searchSite(const std::string& siteName);
	RunPtr getRun(UID id);
	RunPtr searchRun(const std::string& siteName, const UID& runNumber);
	// @}

	// @{ Primary key generation
	virtual UID allocIdRange(const char* tableName, size_t count);
	UID allocIdRangeForConfigs(size_t count);
	UID allocIdRangeForStructureVersions(size_t count);
	UID allocIdRangeForVersions(size_t count);
	UID allocIdRangeForClasses(size_t count);
	UID allocIdRangeForClassInfos(size_t count);
	UID allocIdRangeForAttributes(size_t count);
	UID allocIdRangeForSites(size_t count);
	UID allocIdRangeForRuns(size_t count);
	UID allocIdRangeForRunConditions(size_t count);
	virtual void allocIdRanges(const UID (&counts)[TableCount], UID (&ids)[TableCount]);
	// @}

	// @{ individual entity persistence
	void insert(Config& config);
	void insert(ConfigStructureVersion& structVersion);
	void insert(ConfigVersion& version);
	void insert(Site& site);
	void insert(Run& run);
	void insert(RunCondition& condition);
	// @}

	// @{ entity batch persistence
	void insert(ConfigBatch& batch);
	void insert(ConfigStructureVersionBatch& batch);
	void insert(ConfigVersionBatch& batch);
	void insert(ConfigClassBatch& batch);
	void insert(ConfigClassInfoBatch& batch);
	void insert(ConfigAttributeBatch& batch);
	void insert(SiteBatch& batch);
	// @}

	void insertFullConfig(ConfigStructureVersionBatch& sversions,
			ConfigVersionBatch& versions,
			ConfigClassBatch& classes,
			ConfigClassInfoBatch& classInfos,
			ConfigAttributeBatch& attributes);

	// @{ removal methods
	void remove(Config& config);
	void remove(ConfigStructureVersion& structVersion);
	void remove(ConfigVersion& version);
	void remove(Site& site);
	void remove(Run& run);
	void remove(RunCondition& condition);
	// @}

	// @{ update methods
	void rename(Config& config, const std::string& newName);
	void rename(ConfigVersion& version, const std::string& newName);
	void rename(Site& site, const std::string & newName);
	// @}

protected:
	virtual std::auto_ptr<TableBuilder> createTableBuilder() = 0;
	void createTable(TableBuilder& builder);
	std::string formatPrivileges(int privilege);

	//@{ Single entity construction from a Statement
	ConfigDatabasePtr loadDatabase(SOCI::Statement& stmt);
	ConfigPtr loadConfig(SOCI::Statement& stmt);
	ConfigStructureVersionPtr loadStructure(SOCI::Statement& stmt);
	ConfigVersionPtr loadVersion(SOCI::Statement& stmt);
	SitePtr loadSite(SOCI::Statement& stmt);
	RunPtr loadRun(SOCI::Statement& stmt);
	//@}

	//@{ Multiple entity construction from a Statement
	ConfigPtrSetPtr loadConfigs(SOCI::Statement& stmt);
	ConfigStructureVersionPtrSetPtr loadStructures(SOCI::Statement& stmt);
	ConfigVersionPtrSetPtr loadVersions(SOCI::Statement& stmt);
	ConfigClassPtrSetPtr loadClasses(SOCI::Statement& stmt);
	ConfigClassInfoPtrSetPtr loadClassInfos(SOCI::Statement& stmt);
	ConfigAttributePtrSetPtr loadAttributes(SOCI::Statement& stmt);
	SitePtrSetPtr loadSites(SOCI::Statement& stmt);
	RunPtrSetPtr loadRuns(SOCI::Statement& stmt);
	RunConditionPtrSetPtr loadRunConditions(SOCI::Statement& stmt);
	//@}

private:
	ConfigDao(const ConfigDao& rvalue);
	ConfigDao& operator=(const ConfigDao& rvalue);
};
}
}
}

#endif // CCFG_IO_DB_CONFIGDAO_H_INCLUDED
