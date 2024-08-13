#include "TableBuilderTest.h"
#include "CCfg/Io/Db/MySqlTableBuilder.h"
#include "tut.h"

using namespace CCfg::Io::Db;

namespace tut
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct MySqlTableBuilderTestData
{
	MySqlTableBuilderTestData()
	{
	}
	~MySqlTableBuilderTestData()
	{
	}
};

typedef test_group<MySqlTableBuilderTestData> TestGroup;
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
	MySqlTableBuilder builder;
	std::string query = createQuery1(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tid INTEGER UNSIGNED AUTO_INCREMENT NOT NULL,\n"
		"\tCONSTRAINT Pk_Test\n"
		"\t\tPRIMARY KEY(id)\n"
		") ENGINE=InnoDB;");
}

template<> template<>
void TestObject::test<2>()
{
	MySqlTableBuilder builder;
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
		") ENGINE=InnoDB;"
	);
}

template<> template<>
void TestObject::test<3>()
{
	MySqlTableBuilder builder;
	std::string query = createQuery3(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tname VARCHAR(127) NOT NULL,\n"
		"\tCONSTRAINT Uk_Test\n"
		"\t\tUNIQUE(name)\n"
		") ENGINE=InnoDB;");
}

template<> template<>
void TestObject::test<4>()
{
	MySqlTableBuilder builder;
	std::string query = createQuery4(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tx DECIMAL(15,14),\n"
		"\ta REAL NOT NULL,\n"
		"\tCONSTRAINT Ck_Test\n"
		"\t\tCHECK(a BETWEEN 0 AND 1e-10)\n"
		") ENGINE=InnoDB;");
}

///////////////////////////////////////////////////////////////////////////////
//Test Group instanciation
///////////////////////////////////////////////////////////////////////////////

TestGroup tgMySqlTableBuilder("CCfg::Io::Db::MySqlTableBuilder");

}
