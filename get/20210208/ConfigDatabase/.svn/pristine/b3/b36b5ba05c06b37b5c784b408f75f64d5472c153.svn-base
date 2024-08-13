#include "MySqlTableBuilder.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io:Db::MySqlTableBuilder class methods body
////////////////////////////////////////////////////////////////////////////////

const char* const MySqlTableBuilder::TYPE_NAMES[COLUMN_TYPES_COUNT] =
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
const char* const MySqlTableBuilder::REFERENCE_OPTION_NAMES[REFERENCE_OPTIONS_COUNT] =
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
MySqlTableBuilder::MySqlTableBuilder()
{
}

/**
 * Destructor.
 */
MySqlTableBuilder::~MySqlTableBuilder()
{
}

/*
 * (see inherited doc)
 */
const char* MySqlTableBuilder::getTypeName(ColumnType type) const
{
	return TYPE_NAMES[type];
}

/*
 * (see inherited doc)
 */
const char* MySqlTableBuilder::getReferenceOptionName(ReferenceOption option) const
{
	return REFERENCE_OPTION_NAMES[option];
}

/*
 * (see inherited doc)
 */
const char* MySqlTableBuilder::getCurrentDateTimeSymbol() const
{
	return "CURRENT_TIMESTAMP";
}

/*
 * (see inherited doc)
 */
void MySqlTableBuilder::afterBuild(std::ostringstream& query)
{
	query << " ENGINE=InnoDB";
}

}
}
}
