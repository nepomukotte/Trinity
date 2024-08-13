#ifndef CCFG_IO_DB_MYSQLTABLEBUILDER_H_INCLUDED
#define CCFG_IO_DB_MYSQLTABLEBUILDER_H_INCLUDED

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
 * This subclass defines constants and methods needed to generate MySQL compliant
 * "CREATE TABLE" queries.
 *
 * @author fchateau
 */
class MySqlTableBuilder : public TableBuilder
{
	static const char* const TYPE_NAMES[COLUMN_TYPES_COUNT];
	static const char* const REFERENCE_OPTION_NAMES[REFERENCE_OPTIONS_COUNT];
public:
	MySqlTableBuilder();
	virtual ~MySqlTableBuilder();

protected:
	const char* getTypeName(ColumnType type) const;
	const char* getReferenceOptionName(ReferenceOption option) const;
	const char* getCurrentDateTimeSymbol() const;
	void afterBuild(std::ostringstream& query);

private:
	MySqlTableBuilder(const MySqlTableBuilder& rvalue);
	MySqlTableBuilder& operator=(const MySqlTableBuilder& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_MYSQLTABLEBUILDER_H_INCLUDED
