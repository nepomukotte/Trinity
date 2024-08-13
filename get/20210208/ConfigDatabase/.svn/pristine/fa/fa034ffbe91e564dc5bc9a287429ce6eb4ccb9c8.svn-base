#include "TableBuilderTest.h"
#include "CCfg/Io/Db/OracleTableBuilder.h"
#include "tut.h"

using namespace CCfg::Io::Db;

namespace tut
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct OracleTableBuilderTestData
{
	OracleTableBuilderTestData()
	{
	}
	~OracleTableBuilderTestData()
	{
	}
};

typedef test_group<OracleTableBuilderTestData> TestGroup;
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
	OracleTableBuilder builder;
	std::string query = createQuery1(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tid NUMBER(10) NOT NULL,\n"
		"\tCONSTRAINT Pk_Test\n"
		"\t\tPRIMARY KEY(id)\n"
		");\n"
		"CREATE SEQUENCE db.Test_id;\n"
		"CREATE TRIGGER db.Test_id_Gen\n"
		"BEFORE INSERT ON db.Test\n"
		"FOR EACH ROW\n"
		"BEGIN\n"
		"\tIF :new.id IS NOT NULL THEN\n"
		"\t\tRAISE_APPLICATION_ERROR(-20000, 'ID cannot be specified');\n"
		"\tELSE\n"
		"\t\tSELECT db.Test_id.NEXTVAL\n"
		"\t\tINTO :new.id\n"
		"\t\tFROM dual;\n"
		"\tEND IF;\n"
		"END;"
	);
}

template<> template<>
void TestObject::test<2>()
{
	OracleTableBuilder builder;
	std::string query = createQuery2(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tdt DATE DEFAULT SYSDATE NOT NULL,\n"
		"\tfk NUMBER(10),\n"
		"\tCONSTRAINT Pk_Test\n"
		"\t\tPRIMARY KEY(dt),\n"
		"\tCONSTRAINT Fk_Test_Other\n"
		"\t\tFOREIGN KEY(fk)\n"
		"\t\tREFERENCES db.Other(id)\n"
		"\t\tON DELETE SET NULL\n"
		");"
	);
}

template<> template<>
void TestObject::test<3>()
{
	OracleTableBuilder builder;
	std::string query = createQuery3(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tname VARCHAR2(127) NOT NULL,\n"
		"\tCONSTRAINT Uk_Test\n"
		"\t\tUNIQUE(name)\n"
		");");
}

template<> template<>
void TestObject::test<4>()
{
	OracleTableBuilder builder;
	std::string query = createQuery4(builder);
	ensure_equals(query,
		"CREATE TABLE db.Test(\n"
		"\tx NUMBER(15,14),\n"
		"\ta FLOAT NOT NULL,\n"
		"\tCONSTRAINT Ck_Test\n"
		"\t\tCHECK(a BETWEEN 0 AND 1e-10)\n"
		");");
}


///////////////////////////////////////////////////////////////////////////////
//Test Group instanciation
///////////////////////////////////////////////////////////////////////////////

TestGroup tgOracleTableBuilder("CCfg::Io::Db::OracleTableBuilder");

}
