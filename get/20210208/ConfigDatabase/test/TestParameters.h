#ifndef TESTPARAMETERS_H_
#define TESTPARAMETERS_H_

#include <string>

/**
 * This class parses command line arguments to find connection parameters for
 * each database.
 * The command line must have 3 strings:
 * - the first is the oracle connection string.
 * - the second is the mysql connection string.
 * - the third is the postgresql connection string.
 * - the fourth is the SQLite3 connection string.
 * The connection strings are using the same conventions than the SOCI library,
 * except that the Oracle service must not use a TNS service but declare a
 * service on the fly, like this: "service=host:port/sid"
 */
class TestParameters
{
	static TestParameters parameters;

	TestParameters() {}
	~TestParameters() {}

public:
	std::string connectionString[4];
	std::string host[4];
	std::string port[4];
	std::string user[4];
	std::string passwd[4];
	std::string sid[4];
	std::string schema[4];

	static TestParameters& instance() { return parameters; }
	static void init(const std::string& oracleConnString, const std::string& mysqlConnString,
			const std::string& postgresqlConnString, const std::string & sqlite3ConnString);
};

#endif /*TESTPARAMETERS_H_*/
