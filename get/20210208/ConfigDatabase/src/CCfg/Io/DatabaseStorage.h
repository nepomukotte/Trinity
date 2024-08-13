#ifndef CCFG_IO_DATABASESTORAGE_H_INCLUDED
#define CCFG_IO_DATABASESTORAGE_H_INCLUDED

#include "CCfg/Io/Storage.h"
#include "Db/ConfigDao.h"
#include <boost/shared_ptr.hpp>
#include <string>
#include <memory>
#include <ctime>

namespace CCfg
{
namespace Io
{
// forward declarations
class AbstractLoginCallback;
class Url;

typedef boost::shared_ptr<AbstractLoginCallback> LoginCallbackPtr;

enum DbmsType {Oracle = 0, MySql, PostgreSql, Sqlite3, Invalid = 0xFFFFFFFF};

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is the implementation of the Storage interface for databases. It
 * represents a database connection and an associated current configuration version.
 * The associated configuration is optional, one can use this class only as a
 * "ConfigDao factory".
 *
 * @author fchateau
 */
class DatabaseStorage : public Storage
{
public:
	// classes
	class ConnectionError;

	// static properties
	static const int DEFAULT_PORT[4];
	static const char* const SCHEME[4];
	static const char VERSION_DATE_KEY[];

	static DbmsType dbmsFromString(const std::string& dbmsStr);
	static const char* dbmsToString(DbmsType dbmsType);

	static void fillUrl(Url& url, DbmsType dbms,
			const std::string& user, const std::string& passwd,
			const std::string& host, const std::string& port,
			const std::string& sid,
			const std::string& schema = std::string(),
			const std::string& config = std::string(),
			const std::tm* date = 0);

	static std::string formatUrl(DbmsType dbms,
			const std::string& user, const std::string& passwd,
			const std::string& host, const std::string& port,
			const std::string& sid,
			const std::string& schema = std::string(),
			const std::string& config = std::string(),
			const std::tm* date = 0);

	// @{ construction
	DatabaseStorage();
	DatabaseStorage(DbmsType dbmsType,
			const std::string& host, const std::string& port,
			const std::string& user, const std::string& pwd,
			const std::string& sid, const std::string& schema);
	~DatabaseStorage();
	// @}
	// @{ opening
	void open();
	void close();
	bool isOpened() const;
	// @}
	// @{ accessors
	const std::string& getUrl() const;

	Db::ConfigDao& getDao() { return *dao.get(); }
	const LoginCallbackPtr& getLoginCallback() const { return loginCallback; }
	DbmsType getDbmsType() const { return dbmsType; }
	const std::string& getHost() const { return host; }
	const std::string& getPort() const { return port; }
	const std::string& getUser() const { return user; }
	const std::string& getPassword() const { return pwd; }
	const std::string& getSID() const { return sid; }
	const std::string& getSchema() const { return schema; }
	const std::string& getCurrentConfigName() const { return configName; }
	const std::tm& getCurrentConfigDate() const { return configDate; }

	void setLoginCallback(const LoginCallbackPtr& callback) { loginCallback = callback; }
	void setDbmsType(DbmsType newDbms) { dbmsType = newDbms; }
	void setHost(const std::string& newHost) { host = newHost; }
	void setPort(const std::string& newPort) { port = newPort; }
	void setUser(const std::string& newUser) { user = newUser; }
	void setPassword(const std::string& newPwd) { pwd = newPwd; }
	void setSID(const std::string& newSID) { sid = newSID; }
	void setSchema(const std::string& newSchema) { schema = newSchema; }
	void setCurrentConfigName(const std::string& curConfigName) { configName = curConfigName; }
	void setCurrentConfigDate(const std::tm& curConfigDate) { configDate = curConfigDate; }
	// @}

private:
	DatabaseStorage(const DatabaseStorage& rvalue);
	DatabaseStorage& operator=(const DatabaseStorage& rvalue);

	void checkLoginInfos();

	// @{
	std::auto_ptr<Db::ConfigDao> dao;
	LoginCallbackPtr loginCallback;
	// @}
	// @{ storage properties
	DbmsType dbmsType;
	std::string host;
	std::string port;
	std::string user;
	std::string pwd;
	std::string sid;
	std::string schema;
	std::string configName;
	std::tm configDate;
	mutable std::string url;
	// @}
};

/**
 * This exception signals an error while connecting to the database storage.
 */
class DatabaseStorage::ConnectionError : public Storage::OpenError
{
public:
	ConnectionError(const std::string& message);
	ConnectionError(const ConnectionError& rvalue);
	~ConnectionError() throw() {}
	ConnectionError& operator=(const ConnectionError& rvalue);

private:
	std::string formatMessage() const;
};

}
}

#endif // CCFG_IO_DATABASESTORAGE_H_INCLUDED
