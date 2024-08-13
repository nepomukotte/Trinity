#include "TestParameters.h"
#include "tut.h"
#include "tut_reporter.h"
#include <iostream>
#include <fstream>

namespace tut
{
	test_runner_singleton runner;
}

using tut::runner;
using tut::reporter;

int main(int /* argc */, char* /* argv */[])
{
	std::ifstream configFile("connections.txt");
	if (configFile.fail())
	{
		std::cout << "\n\nThe program requires a text file named \"connections.txt\" "
			"containing 4 lines, in this order: \n"
			"1/ the oracle connection string, using the following format:\n"
			"\tservice=<host>:<port>/<sid> user=<user> password=<passwd>\n"
			"2/ the mysql connection string:\n"
			"\thost=<host> port=<port> user=<user> password=<passwd> db=<db>\n"
			"3/ the postgresql connection string:\n"
			"\thost=<host> port=<port> user=<user> password=<passwd> db=<db>\n"
			"4/ the SQLite3 connection string:\n"
			"\t<path>\n"
			"The port is optional if using the default port.\n";
		return -1;
	}
	std::string oracleConnection, mysqlConnection, pgsqlConnection, sqlite3Connection;
	getline(configFile, oracleConnection);
	getline(configFile, mysqlConnection);
	getline(configFile, pgsqlConnection);
	getline(configFile, sqlite3Connection);
	TestParameters::init(oracleConnection, mysqlConnection, pgsqlConnection, sqlite3Connection);

	// instantiate and set the callback that displays test run reports
	reporter callback; // use std::cout by default to output report.
	runner.get().set_callback(&callback);

	// run all tests in all groups
	runner.get().run_tests();

	return 0;
}
