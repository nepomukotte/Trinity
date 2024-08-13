#ifndef CCFG_IO_DB_ORACLETABLEBUILDER_H_INCLUDED
#define CCFG_IO_DB_ORACLETABLEBUILDER_H_INCLUDED

#include "TableBuilder.h"
#include <vector>
#include <string>

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
 * This subclass defines constants and methods needed to generate Oracle compliant
 * "CREATE TABLE" queries.
 *
 * @author fchateau
 */
class OracleTableBuilder : public TableBuilder
{
	static const char* const TYPE_NAMES[COLUMN_TYPES_COUNT];
	static const char* const REFERENCE_OPTION_NAMES[REFERENCE_OPTIONS_COUNT];
	std::vector<std::string> autoIncrementColumns;
	std::vector<std::string> sideQueries;

public:
	OracleTableBuilder();
	virtual ~OracleTableBuilder();
	std::vector<std::string> getSideQueries();

protected:
	const char* getTypeName(ColumnType type) const;
	const char* getReferenceOptionName(ReferenceOption option) const;
	const char* getCurrentDateTimeSymbol() const;
	void beforeBuild(std::ostringstream& query);
	void afterBuild(std::ostringstream& query);

private:
	OracleTableBuilder(const OracleTableBuilder& rvalue);
	OracleTableBuilder& operator=(const OracleTableBuilder& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_ORACLETABLEBUILDER_H_INCLUDED
