#include "DatabaseStorage.h"
#include "AbstractLoginCallback.h"
#include "CCfg/Io/Db/MySqlConfigDao.h"
#ifdef CCFG_USE_ORACLE
#include "CCfg/Io/Db/OracleConfigDao.h"
#endif
#ifdef CCFG_USE_POSTGRESQL
#include "CCfg/Io/Db/PostgreSqlConfigDao.h"
#endif
#ifdef CCFG_USE_SQLITE3
#include "CCfg/Io/Db/Sqlite3ConfigDao.h"
#endif
#include "CCfg/Io/Url.h"
#include "Utils/Date/StdTmHelpers.h"
#include <soci.h>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <vector>
#include <algorithm>

namespace CCfg
{
namespace Io
{
using namespace Db;
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::DatabaseStorage class methods body
////////////////////////////////////////////////////////////////////////////////

const int DatabaseStorage::DEFAULT_PORT[4] = {1521, 3306, 5432, -1};
const char* const DatabaseStorage::SCHEME[4] = {"oracle", "mysql", "postgresql", "sqlite3"};
const char DatabaseStorage::VERSION_DATE_KEY[] = "date";

/**
 * Returns the Dbms enum corresponding to a database scheme string.
 */
DbmsType DatabaseStorage::dbmsFromString(const std::string& dbmsStr)
{
	DbmsType dbms;
	if(dbmsStr == DatabaseStorage::SCHEME[Oracle])
	{
		dbms = Oracle;
	}
	else if(dbmsStr == DatabaseStorage::SCHEME[MySql])
	{
		dbms = MySql;
	}
	else if(dbmsStr == DatabaseStorage::SCHEME[PostgreSql])
	{
		dbms = PostgreSql;
	}
	else if(dbmsStr == DatabaseStorage::SCHEME[Sqlite3])
	{
		dbms = Sqlite3;
	}
	else
	{
		dbms = Invalid;
	}
	return dbms;
}

/**
 * Returns the scheme corresponding to a Dbms enum value.
 */
const char* DatabaseStorage::dbmsToString(DbmsType dbmsType)
{
	if(dbmsType != Oracle and dbmsType != MySql and dbmsType != PostgreSql and dbmsType != Sqlite3)
	{
		return "invalid";
	}
	else
	{
		return DatabaseStorage::SCHEME[dbmsType];
	}
}

/**
 * Initializes an URL with the informations composing a database URL.
 */
void DatabaseStorage::fillUrl(Url& url, DbmsType dbms,
		const std::string& user, const std::string& passwd,
		const std::string& host, const std::string& port,
		const std::string& sid, const std::string& schema,
		const std::string& config, const std::tm* date)
{
	url.setScheme(dbmsToString(dbms));
	url.setUser(user);
	if(not passwd.empty())	url.setPassword(passwd);
	url.setHost(host);
	if(not port.empty())	url.setPort(boost::lexical_cast<int>(port));
	if(not sid.empty())		url.appendPath(sid);
	if(not schema.empty())	url.appendPath(schema);
	if(not config.empty())	url.appendPath(config);
	if(date != 0)
	{
		std::ostringstream str;
		str << *date;
		url.addQueryItem(DatabaseStorage::VERSION_DATE_KEY, str.str());
	}
}

/**
 * Creates a string URL from the informations composing a database URL.
 */
std::string DatabaseStorage::formatUrl(DbmsType dbms,
		const std::string& user, const std::string& passwd,
		const std::string& host, const std::string& port,
		const std::string& sid, const std::string& schema,
		const std::string& config, const std::tm* date)
{
	Url url;
	DatabaseStorage::fillUrl(url, dbms, user, passwd, host, port, sid, schema, config, date);
	return url.format();
}

/**
 * Default constructor.
 */
DatabaseStorage::DatabaseStorage()
	: dbmsType(Oracle), host(), port(), user(), pwd(),
	sid(), schema(), configName(), configDate(), url()
{
}

/**
 * Full constructor.
 */
DatabaseStorage::DatabaseStorage(DbmsType dbmsType,
	const std::string& host, const std::string& port,
	const std::string& user, const std::string& pwd,
	const std::string& sid, const std::string& schema)
	: dbmsType(dbmsType), host(host), port(port), user(user), pwd(pwd),
	sid(sid), schema(schema), configName(), configDate(), url()
{
	assert(dbmsType == Oracle or dbmsType == MySql or dbmsType == PostgreSql or dbmsType == Sqlite3);
}

/**
 * Destructor.
 */
DatabaseStorage::~DatabaseStorage()
{
}

/*
 * (see inherited doc)
 */
void DatabaseStorage::open()
{
	if (host.empty() and dbmsType != Sqlite3)
	{
		throw ConnectionError("host not specified.");
	}

	std::vector< DbmsType > enabledTypes;
	enabledTypes.push_back(MySql);
#if defined(CCFG_USE_ORACLE)
	enabledTypes.push_back(Oracle);
#endif
#if defined(CCFG_USE_POSTGRESQL)
	enabledTypes.push_back(PostgreSql);
#endif
#if defined(CCFG_USE_SQLITE3)
	enabledTypes.push_back(Sqlite3);
#endif
	assert(std::find(enabledTypes.begin(), enabledTypes.end(), dbmsType) != enabledTypes.end());
	try
	{
		if ((user.empty() or pwd.empty()) and dbmsType != Sqlite3)
		{
			if(not loginCallback->requestLogin(dbmsType, host, port, sid, user, pwd))
			{
				return;
			}
		}
		if(dbmsType == MySql)
		{
			MySqlConfigDao* mysqlDao = new MySqlConfigDao;
			dao.reset(mysqlDao);
			mysqlDao->connect(host, port, user, pwd, schema);
		}
#ifdef CCFG_USE_ORACLE
		else if(dbmsType == Oracle)
		{
			OracleConfigDao* oracleDao = new OracleConfigDao;
			dao.reset(oracleDao);
			oracleDao->connect(host, port, sid, user, pwd);
			if(not schema.empty())
			{
				oracleDao->setWorkingDatabase(schema);
			}
			else
			{
				schema = oracleDao->getWorkingDatabase();
			}
		}
#endif
#ifdef CCFG_USE_POSTGRESQL
		else if (dbmsType == PostgreSql)
		{
			PostgreSqlConfigDao* postgresqlDao = new PostgreSqlConfigDao;
			dao.reset(postgresqlDao);
			postgresqlDao->connect(host, port, sid, user, pwd);
			if(not schema.empty())
			{
				postgresqlDao->setWorkingDatabase(schema);
			}
			else
			{
				schema = postgresqlDao->getWorkingDatabase();
			}
		}
#endif
#ifdef CCFG_USE_SQLITE3
		else if (dbmsType == Sqlite3)
		{
			Sqlite3ConfigDao* sqlite3Dao = new Sqlite3ConfigDao;
			dao.reset(sqlite3Dao);
			sqlite3Dao->connect(sid);
			if(not schema.empty())
			{
				sqlite3Dao->setWorkingDatabase(schema);
			}
			else
			{
				schema = sqlite3Dao->getWorkingDatabase();
			}
		}
#endif
	}
	catch(SOCI::SOCIError& e)
	{
		pwd.clear();
		loginCallback->onLoginFailed(dbmsType, host, port, sid, user);
		throw ConnectionError(e.what());
	}
}

/*
 * (see inherited doc)
 */
void DatabaseStorage::close()
{
	dao.reset();
}

/*
 * (see inherited doc)
 */
bool DatabaseStorage::isOpened() const
{
	return dao.get() != 0 and dao->isConnected();
}

/*
 * (see inherited doc)
 */
const std::string& DatabaseStorage::getUrl() const
{
	url = formatUrl(dbmsType, user, pwd, host, port, sid, schema, configName,
			(configDate.tm_year != 0) ? &configDate : 0);
	return url;
}

////////////////////////////////////////////////////////////////////////////////
// classname class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor from an error message.
 */
DatabaseStorage::ConnectionError::ConnectionError(const std::string& message) : OpenError(message)
{
}

/**
 * Copy constructor.
 */
DatabaseStorage::ConnectionError::ConnectionError(const ConnectionError& rvalue)
	: OpenError(rvalue)
{
}

/**
 * Assignment operator.
 */
DatabaseStorage::ConnectionError&
DatabaseStorage::ConnectionError::operator=(const ConnectionError& rvalue)
{
	OpenError::operator=(rvalue);
	return *this;
}

/*
 * (see inherited doc)
 */
std::string DatabaseStorage::ConnectionError::formatMessage() const
{
	std::ostringstream errMsg;
	errMsg << "Error connecting to the database: " << getMessage();
	return errMsg.str();
}

}
}
