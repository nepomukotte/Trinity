#include "TableBuilderTest.h"
#include "CCfg/Io/Db/Sqlite3TableBuilder.h"
#include "tut.h"

using namespace CCfg::Io::Db;

namespace tut
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct Sqlite3TableBuilderTestData
{
	Sqlite3TableBuilderTestData()
	{
	}
	~Sqlite3TableBuilderTestData()
	{
	}
};

typedef test_group<Sqlite3TableBuilderTestData> TestGroup;
typedef TestGroup::object TestObject;


///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 *
 */
template<> template<>
void TestObject::test<1>()
{
	Sqlite3TableBuilder builder;
	std::string query = createQuery1(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tid INTEGER UNSIGNED AUTO_INCREMENT NOT NULL,\n"
		"\tCONSTRAINT Pk_Test\n"
		"\t\tPRIMARY KEY(id)\n"
		");");
}

template<> template<>
void TestObject::test<2>()
{
	Sqlite3TableBuilder builder;
	std::string query = createQuery2(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tdt TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL,\n"
		"\tfk INTEGER,\n"
		"\tCONSTRAINT Pk_Test\n"
		"\t\tPRIMARY KEY(dt),\n"
		"\tCONSTRAINT Fk_Test_Other\n"
		"\t\tFOREIGN KEY(fk)\n"
		"\t\tREFERENCES Other(id)\n"
		"\t\tON DELETE RESTRICT\n"
		"\t\tON UPDATE CASCADE\n"
		");"
	);
}

template<> template<>
void TestObject::test<3>()
{
	Sqlite3TableBuilder builder;
	std::string query = createQuery3(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tname VARCHAR(127) NOT NULL,\n"
		"\tCONSTRAINT Uk_Test\n"
		"\t\tUNIQUE(name)\n"
		");");
}

template<> template<>
void TestObject::test<4>()
{
	Sqlite3TableBuilder builder;
	std::string query = createQuery4(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tx DECIMAL(15,14),\n"
		"\ta REAL NOT NULL,\n"
		"\tCONSTRAINT Ck_Test\n"
		"\t\tCHECK(a BETWEEN 0 AND 1e-10)\n"
		");");
}

///////////////////////////////////////////////////////////////////////////////
//Test Group instanciation
///////////////////////////////////////////////////////////////////////////////

TestGroup tgSqlite3TableBuilder("CCfg::Io::Db::Sqlite3TableBuilder");

}
