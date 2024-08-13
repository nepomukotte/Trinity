#include "TestParameters.h"
#include "CCfg/Io/DatabaseStorage.h"
#include <boost/tokenizer.hpp>
#include <string>

using CCfg::Io::Oracle;
using CCfg::Io::MySql;
using CCfg::Io::PostgreSql;
using CCfg::Io::Sqlite3;

TestParameters TestParameters::parameters;

void TestParameters::init(const std::string& oracleConnString, const std::string& mysqlConnString,
		const std::string& pgsqlConnString, const std::string & sqlite3ConnString)
{
	parameters.connectionString[Oracle] = oracleConnString;
	parameters.connectionString[MySql] = mysqlConnString;
	parameters.connectionString[PostgreSql] = pgsqlConnString;
	parameters.connectionString[Sqlite3] = sqlite3ConnString;

	typedef boost::char_separator<char> Separator;
	typedef boost::tokenizer<Separator> Tokenizer;

	Separator sep(" \t");
	for(size_t c = 0; c < 4; ++c)
	{
		Tokenizer tok(parameters.connectionString[c], sep);
		for(Tokenizer::iterator cur = tok.begin(), end = tok.end(); cur != end; ++cur)
		{
			const std::string& key_value = *cur;
			size_t equalSignIndex = key_value.find_first_of('=');
			std::string key = key_value.substr(0, equalSignIndex);
			std::string value = key_value.substr(equalSignIndex + 1, std::string::npos);
			if(key == "host")
			{
				parameters.host[c] = value;
			}
			else if(key == "port")
			{
				parameters.port[c] = value;
			}
			else if(key == "user")
			{
				parameters.user[c] = value;
			}
			else if(key == "password")
			{
				parameters.passwd[c] = value;
			}
			else if(key == "db")
			{
				parameters.schema[c] = value;
			}
			else if (key =="dbname")
			{
				parameters.sid[c] = value;
			}
			else if (key =="schema")
			{
				parameters.schema[c] = value;
			}
			else if(key == "service")
			{
				size_t colonIndex = value.find_first_of(':');
				size_t slashIndex = value.find_first_of('/');
				parameters.host[c] = value.substr(0, std::min(colonIndex, slashIndex));
				if(colonIndex != std::string::npos)
				{
					parameters.port[c] =  value.substr(colonIndex + 1, slashIndex - (colonIndex + 1));
				}
				if(slashIndex != std::string::npos)
				{
					parameters.sid[c] = value.substr(slashIndex + 1, std::string::npos);
				}
			}
		}
	}
	if(parameters.port[Oracle].empty())
	{
		parameters.port[Oracle] = "1521";
	}
	if(parameters.port[MySql].empty())
	{
		parameters.port[MySql] = "3306";
	}
	if(parameters.port[PostgreSql].empty())
	{
		parameters.port[PostgreSql] = "5432";
	}
	if(parameters.schema[Oracle].empty())
	{
		parameters.schema[Oracle] = parameters.user[Oracle];
	}
}
