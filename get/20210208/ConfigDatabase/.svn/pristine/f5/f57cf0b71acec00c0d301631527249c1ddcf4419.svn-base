#ifndef CCFG_IO_DB_TABLEBUILDER_H_INCLUDED
#define CCFG_IO_DB_TABLEBUILDER_H_INCLUDED

#include <vector>
#include <sstream>
#include <string>
#include <memory>

namespace CCfg
{
namespace Io
{
namespace Db
{
// enums
enum ColumnType {Integer = 0, ///< Integer represented with 32 bits at least
	Decimal,	///< Fixed point number. Use length and decimals to describe the number of digits.
	Real,		///< Floating point number, with exponent range in [-125;125] at least.
	Boolean,	///< Boolean: 0 or 1.
	Varchar,	///< Variable length character string. The maximum length must be specified in the length parameter.
	Char,		///< Fixed length character string. The length must be specified in the length parameter.
	DateTime,	///< Date and time with 1 second precision.
	Date,		///< Date with 1 day precision.
	Time		///< Time (without date) with 1 second precision.
};
enum ReferenceOption {Restrict = 0,	///< forbids the action if the record is still referenced.
	Cascade,	///< for update actions it updates the foreign key of referencing records. For drop actions it also drops all referencing records.
	SetNull,	///< sets the foreign key of referencing records to Null.
	NoAction,	///< do nothing
	SetDefault	///< resets the foreign key of referencing records to its default value, if any.
};
enum ColumnSign {SignNA,	///< the concept of sign is not application to the specified type.
	Signed,		///< the type has a sign: positive and negative values can be represented.
	Unsigned	///< the type is unsigned: only positive values can be represented.
};
enum NullPolicy {Null,	///< Null values are accepted
	NotNull	///< Null values are forbidden
};
enum InsertPolicy {Manual,	///< the value must be manually specified.
	AutoIncrement,	///< the value will be automatically generated using an incremental integer. Manual definition may be restricted.
	CurrentDateTime	///< the value will be automatically generated using the current date.
};
#define COLUMN_TYPES_COUNT	Time+1
#define REFERENCE_OPTIONS_COUNT		SetDefault+1
#define UNDEFINED	0xFFFFFFFFu

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class is useful to generate "CREATE TABLE" queries easily and in a RDMS
 * independent way.
 * To use this class, first the table must be described with methods like
 * setTableName, addColumn, etc... and then the query can be build with
 * buildQuery.
 * The order in which the description methods are called does not matter, but in
 * any case, setTableName and setDatabase must be called.
 *
 * @author fchateau
 */
class TableBuilder
{
protected:
	// inner classes;
	struct ColumnDesc;
	struct ConstraintDesc;
	struct ForeignKeyDesc;
	struct UniqueDesc;
	struct CheckDesc;

	// attributes
	const char* dbName;
	const char* tableName;
	std::vector<ColumnDesc> columns;
	const char* primaryKeyColumns;
	const char* pkConstraintName;
	std::vector<ForeignKeyDesc> foreignKeys;
	std::vector<UniqueDesc> uniqueColumns;
	std::vector<CheckDesc> checks;

public:
	TableBuilder();
	virtual ~TableBuilder();

	void setDatabase(const char* name);
	void setTableName(const char* name);
	void addColumn(const char* name, ColumnType type, ColumnSign sign, NullPolicy nullPolicy, InsertPolicy insertPolicy = Manual);
	void addColumn(const char* name, ColumnType type, ColumnSign sign, unsigned length, unsigned decimals, NullPolicy nullPolicy, InsertPolicy insertPolicy = Manual);
	void addPrimaryKey(const char* constraintName, const char* columns);
	void addForeignKey(const char* constraintName, const char* fkColumns, const char* refTable, const char* refColumns, ReferenceOption onUpdate, ReferenceOption onDelete);
	void addUnique(const char* constraintName, const char* columnsList);
	void addCheck(const char* constraintName, const char* checkExpr);
	std::string buildQuery();
	virtual std::vector<std::string> getSideQueries();

protected:
	virtual std::string formatForeignKey(const ForeignKeyDesc & fkDesc);
	virtual const char* getTypeName(ColumnType type) const = 0;
	virtual const char* getReferenceOptionName(ReferenceOption option) const = 0;
	virtual const char* getCurrentDateTimeSymbol() const = 0;
	virtual void beforeBuild(std::ostringstream& query);
	virtual void afterBuild(std::ostringstream& query);

private:
	TableBuilder(const TableBuilder& rvalue);
	TableBuilder& operator=(const TableBuilder& rvalue);
};


struct TableBuilder::ColumnDesc
{
	const char* name;
	ColumnType type;
	unsigned length;
	unsigned decimals;
	ColumnSign sign;
	NullPolicy nullPolicy;
	InsertPolicy insertPolicy;
};

struct TableBuilder::ConstraintDesc
{
	const char* constraintName;
};

struct TableBuilder::ForeignKeyDesc : TableBuilder::ConstraintDesc
{
	const char* fkColumns;
	const char* refTable;
	const char* refColumns;
	ReferenceOption onUpdate;
	ReferenceOption onDelete;
};

struct TableBuilder::UniqueDesc : TableBuilder::ConstraintDesc
{
	const char* columnsList;
};

struct TableBuilder::CheckDesc : TableBuilder::ConstraintDesc
{
	const char* checkExpr;
};

}
}
}

#endif // CCFG_IO_DB_TABLEBUILDER_H_INCLUDED
