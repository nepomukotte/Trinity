#include "TableBuilderTest.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
std::string buildFullQuery(TableBuilder& builder)
{
	std::string query = builder.buildQuery();
	if(*query.rbegin() != ';')
	{
		query += ';';
	}

	const std::vector<std::string>& sideQueries = builder.getSideQueries();
	for(size_t i = 0; i < sideQueries.size(); ++i)
	{
		query += '\n';
		query += sideQueries[i];
		if(*query.rbegin() != ';')
		{
			query += ';';
		}
	}
	return query;
}

std::string createQuery1(TableBuilder& builder)
{
	builder.setDatabase("db");
	builder.setTableName("Test");
	builder.addColumn("id", Integer, Unsigned, NotNull, AutoIncrement);
	builder.addPrimaryKey("Pk_Test", "id");
	return buildFullQuery(builder);
}

std::string createQuery2(TableBuilder& builder)
{
	builder.setDatabase("db");
	builder.setTableName("Test");
	builder.addColumn("dt", DateTime, SignNA, NotNull, CurrentDateTime);
	builder.addColumn("fk", Integer, Signed, Null);
	builder.addPrimaryKey("Pk_Test", "dt");
	builder.addForeignKey("Fk_Test_Other", "fk", "Other", "id", Cascade, Restrict);
	return buildFullQuery(builder);
}

std::string createQuery3(TableBuilder& builder)
{
	builder.setDatabase("db");
	builder.setTableName("Test");
	builder.addColumn("name", Varchar, SignNA, 127, UNDEFINED, NotNull);
	builder.addUnique("Uk_Test", "name");
	return buildFullQuery(builder);
}

std::string createQuery4(TableBuilder& builder)
{
	builder.setDatabase("db");
	builder.setTableName("Test");
	builder.addColumn("x", Decimal, Signed, 15, 14, Null);
	builder.addColumn("a", Real, Signed, NotNull);
	builder.addCheck("Ck_Test", "a BETWEEN 0 AND 1e-10");
	return buildFullQuery(builder);
}

}
}
}
