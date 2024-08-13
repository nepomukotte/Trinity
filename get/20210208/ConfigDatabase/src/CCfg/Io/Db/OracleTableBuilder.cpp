#include "OracleTableBuilder.h"
#include <algorithm>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::OracleTableBuilder class methods body
////////////////////////////////////////////////////////////////////////////////

const char* const OracleTableBuilder::TYPE_NAMES[COLUMN_TYPES_COUNT] =
{
	"NUMBER",
	"NUMBER",
	"FLOAT",
	"NUMBER",
	"VARCHAR2",
	"CHAR",
	"DATE",
	"DATE",
	"DATE"
};
const char* const OracleTableBuilder::REFERENCE_OPTION_NAMES[REFERENCE_OPTIONS_COUNT] =
{
	"",
	"CASCADE",
	"SET NULL",
	"",
	""
};

/**
 * Default constructor.
 */
OracleTableBuilder::OracleTableBuilder()
	: autoIncrementColumns(), sideQueries()
{
}

/**
 * Destructor.
 */
OracleTableBuilder::~OracleTableBuilder()
{
}

/*
 * (see inherited doc)
 */
std::vector<std::string> OracleTableBuilder::getSideQueries()
{
	return sideQueries;
}

/*
 * (see inherited doc)
 */
const char* OracleTableBuilder::getTypeName(ColumnType type) const
{
	return TYPE_NAMES[type];
}

/*
 * (see inherited doc)
 */
const char* OracleTableBuilder::getReferenceOptionName(ReferenceOption option) const
{
	return REFERENCE_OPTION_NAMES[option];
}

/*
 * (see inherited doc)
 */
const char* OracleTableBuilder::getCurrentDateTimeSymbol() const
{
	return "SYSDATE";
}

/*
 * (see inherited doc)
 */
void OracleTableBuilder::beforeBuild(std::ostringstream& query)
{
	for(size_t i = 0; i < columns.size(); ++i)
	{
		ColumnDesc& column = columns[i];
		if(column.type == Integer)
		{
			column.length = std::min(10u, column.length);
		}
		if(column.type == Boolean)
		{
			column.length = 1;
		}
		column.sign = Signed;
		if(column.insertPolicy == AutoIncrement)
		{
			autoIncrementColumns.push_back(column.name);
			column.insertPolicy = Manual;
		}
	}
	for(size_t i = 0; i < foreignKeys.size(); ++i)
	{
		ForeignKeyDesc& fkDesc = foreignKeys[i];
		fkDesc.onUpdate = NoAction; // Oracle doesn't support ON UPDATE clause, this value prevents the builder from adding it.
		if(fkDesc.onDelete != Cascade and fkDesc.onDelete != SetNull)
		{
			fkDesc.onDelete = SetNull;
		}
	}
}

/*
 * (see inherited doc)
 */
void OracleTableBuilder::afterBuild(std::ostringstream& query)
{
	for(size_t i = 0; i < autoIncrementColumns.size(); ++i)
	{
		const std::string& columnName = autoIncrementColumns[i];
		std::ostringstream sequenceQuery;
		sequenceQuery <<
			"CREATE SEQUENCE " << dbName << '.' <<
			tableName << '_' << columnName;
		sideQueries.push_back(sequenceQuery.str());

		std::ostringstream triggerQuery;
		triggerQuery <<
			"CREATE TRIGGER " << dbName << '.' <<
				tableName << '_' << columnName << "_Gen\n"
			"BEFORE INSERT ON " << dbName << '.' << tableName << std::endl <<
			"FOR EACH ROW\n"
			"BEGIN\n"
			"\tIF :new.id IS NOT NULL THEN\n"
			"\t\tRAISE_APPLICATION_ERROR(-20000, 'ID cannot be specified');\n"
			"\tELSE\n"
			"\t\tSELECT " << dbName << '.' << tableName << '_' << columnName << ".NEXTVAL\n"
			"\t\tINTO :new.id\n"
			"\t\tFROM dual;\n"
			"\tEND IF;\n"
			"END;";
		sideQueries.push_back(triggerQuery.str());
	}
}

}
}
}
