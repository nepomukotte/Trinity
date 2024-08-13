#ifndef CCFG_IO_DB_SQLITE3TABLEBUILDER_H_INCLUDED
#define CCFG_IO_DB_SQLITE3TABLEBUILDER_H_INCLUDED

#include "TableBuilder.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This subclass defines constants and methods needed to generate SQLite 3 compliant
 * "CREATE TABLE" queries.
 *
 * @author psizun
 */
class Sqlite3TableBuilder : public TableBuilder
{
	static const char* const TYPE_NAMES[COLUMN_TYPES_COUNT];
	static const char* const REFERENCE_OPTION_NAMES[REFERENCE_OPTIONS_COUNT];
public:
	Sqlite3TableBuilder();
	virtual ~Sqlite3TableBuilder();

protected:
	std::string formatForeignKey(const ForeignKeyDesc & fkDesc);
	const char* getTypeName(ColumnType type) const;
	const char* getReferenceOptionName(ReferenceOption option) const;
	const char* getCurrentDateTimeSymbol() const;
	void afterBuild(std::ostringstream& query);

private:
	Sqlite3TableBuilder(const Sqlite3TableBuilder& rvalue);
	Sqlite3TableBuilder& operator=(const Sqlite3TableBuilder& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_SQLITE3TABLEBUILDER_H_INCLUDED
