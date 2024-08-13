#include "TableBuilderTest.h"
#include "CCfg/Io/Db/PostgreSqlTableBuilder.h"
#include "tut.h"
using namespace CCfg::Io::Db;

namespace tut
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct PostgreSqlTableBuilderTestData
{
	PostgreSqlTableBuilderTestData()
	{
	}
	~PostgreSqlTableBuilderTestData()
	{
	}
};

typedef test_group<PostgreSqlTableBuilderTestData> TestGroup;
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
	PostgreSqlTableBuilder builder;
	std::string query = createQuery1(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tid INTEGER NOT NULL,\n"
		"\tCONSTRAINT Pk_Test\n"
		"\t\tPRIMARY KEY(id)\n"
		");\n"
		"CREATE SEQUENCE db.Test_id OWNED BY db.Test.id;\n"
		"ALTER TABLE ONLY db.Test\n"
		"\tALTER COLUMN id SET DEFAULT nextval(\'db.Test_id\');");
}

template<> template<>
void TestObject::test<2>()
{
	PostgreSqlTableBuilder builder;
	std::string query = createQuery2(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tdt TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL,\n"
		"\tfk INTEGER,\n"
		"\tCONSTRAINT Pk_Test\n"
		"\t\tPRIMARY KEY(dt),\n"
		"\tCONSTRAINT Fk_Test_Other\n"
		"\t\tFOREIGN KEY(fk)\n"
		"\t\tREFERENCES db.Other(id)\n"
		"\t\tON DELETE RESTRICT\n"
		"\t\tON UPDATE CASCADE\n"
		");"
	);
}

template<> template<>
void TestObject::test<3>()
{
	PostgreSqlTableBuilder builder;
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
	PostgreSqlTableBuilder builder;
	std::string query = createQuery4(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tx NUMERIC(15,14),\n"
		"\ta REAL NOT NULL,\n"
		"\tCONSTRAINT Ck_Test\n"
		"\t\tCHECK(a BETWEEN 0 AND 1e-10)\n"
		");");
}

///////////////////////////////////////////////////////////////////////////////
//Test Group instantiation
///////////////////////////////////////////////////////////////////////////////

TestGroup tgPostgreSqlTableBuilder("CCfg::Io::Db::PostgreSqlTableBuilder");

}
