#include "Sqlite3TableBuilder.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io:Db::Sqlite3TableBuilder class methods body
////////////////////////////////////////////////////////////////////////////////

const char* const Sqlite3TableBuilder::TYPE_NAMES[COLUMN_TYPES_COUNT] =
{
	"INTEGER",
	"DECIMAL",
	"REAL",
	"BOOLEAN",
	"VARCHAR",
	"CHAR",
	"TIMESTAMP",
	"DATE",
	"TIME"
};
const char* const Sqlite3TableBuilder::REFERENCE_OPTION_NAMES[REFERENCE_OPTIONS_COUNT] =
{
	"RESTRICT",
	"CASCADE",
	"SET NULL",
	"NO ACTION",
	"SET DEFAULT"
};

/**
 * Default constructor.
 */
Sqlite3TableBuilder::Sqlite3TableBuilder()
{
}

/**
 * Destructor.
 */
Sqlite3TableBuilder::~Sqlite3TableBuilder()
{
}

/**
 * Format foreign key as SQLite 3 does not accept database name ('main')
 *  to be prefixed to reference table name.
 */
std::string Sqlite3TableBuilder::formatForeignKey(const ForeignKeyDesc & fkDesc)
{
	std::ostringstream fkConstraint;
	fkConstraint <<
			"FOREIGN KEY(" << fkDesc.fkColumns << ')' << "\n\t\t" <<
			"REFERENCES " << fkDesc.refTable <<
				'(' << fkDesc.refColumns << ')';
	return fkConstraint.str();
}

/*
 * (see inherited doc)
 */
const char* Sqlite3TableBuilder::getTypeName(ColumnType type) const
{
	return TYPE_NAMES[type];
}

/*
 * (see inherited doc)
 */
const char* Sqlite3TableBuilder::getReferenceOptionName(ReferenceOption option) const
{
	return REFERENCE_OPTION_NAMES[option];
}

/*
 * (see inherited doc)
 */
const char* Sqlite3TableBuilder::getCurrentDateTimeSymbol() const
{
	return "CURRENT_TIMESTAMP";
}

/*
 * (see inherited doc)
 */
void Sqlite3TableBuilder::afterBuild(std::ostringstream& query)
{
	;
}

}
}
}
