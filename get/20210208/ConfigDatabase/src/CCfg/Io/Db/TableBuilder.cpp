#include "TableBuilder.h"
#include <iostream>
namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::TableBuilder class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
TableBuilder::TableBuilder()
	: dbName(0), tableName(0), columns(),
	primaryKeyColumns(0), pkConstraintName(0),
	foreignKeys(), uniqueColumns(), checks()
{
}

/**
 * Destructor.
 */
TableBuilder::~TableBuilder()
{
}

/**
 * Defines the database where the table will be created.
 * For the moment, calling this method is mandatory, even if the table should
 * be created in the current database.
 * this may change in a future version.
 */
void TableBuilder::setDatabase(const char* name)
{
	dbName = name;
}

/**
 * Defines the name of the table to create.
 * Calling this method is mandatory.
 */
void TableBuilder::setTableName(const char* name)
{
	tableName = name;
}

/**
 * Adds a column to create.
 * @param name the name of the column.
 * @param type the type of the column.
 * @param sign indicates if the column data is signed or unsigned.
 * @param nullPolicy indicates whether the column should accept null values or not.
 * @param insertPolicy indicates how the value gets generated when inserting a row.
 */
void TableBuilder::addColumn(const char* name, ColumnType type, ColumnSign sign, NullPolicy nullPolicy, InsertPolicy insertPolicy)
{
	ColumnDesc desc;
	desc.name = name;
	desc.type = type;
	desc.sign = sign;
	desc.length = UNDEFINED;
	desc.decimals = UNDEFINED;
	desc.nullPolicy = nullPolicy;
	desc.insertPolicy = insertPolicy;
	columns.push_back(desc);
}

/**
 * Adds a column to the table definition.
 * @param name the name of the column.
 * @param type the type of the column.
 * @param sign indicates if the column data is signed or unsigned.
 * @param length indicates the length of the field.
 * For integers, and decimal types, it represents the maximum total number of
 * digits (including decimals).
 * For real numbers, it represents the maximum exponent (base 10).
 * For strings, it represents the maximum length.
 * For temporal types, the meaning is database specific.
 * @param decimals indicates the maximum number of digits after the decimal point.
 * @param nullPolicy indicates whether the column should accept null values or not.
 * @param insertPolicy indicates how the value gets generated when inserting a row.
 */
void TableBuilder::addColumn(const char* name, ColumnType type, ColumnSign sign, unsigned length, unsigned decimals, NullPolicy nullPolicy, InsertPolicy insertPolicy)
{
	ColumnDesc desc;
	desc.name = name;
	desc.type = type;
	desc.sign = sign;
	desc.length = length;
	desc.decimals = decimals;
	desc.nullPolicy = nullPolicy;
	desc.insertPolicy = insertPolicy;
	columns.push_back(desc);
}

/**
 * Adds a primary key to the table definition.
 * @param constraintName the name of the constraint object in the database.
 * @param columns the name of the column, or the list of comma-separated columns
 * name composing the primary key.
 */
void TableBuilder::addPrimaryKey(const char* constraintName, const char* columns)
{
	pkConstraintName = constraintName;
	primaryKeyColumns = columns;
}

/**
 * Adds a foreign key to the table definition.
 * @param constraintName the name for the constraint.
 * @param fkColumns the name of the column, or the list of comma-separated columns
 * name composing the foreign key.
 * @param refTable the name of the table referenced by the key.
 * @param refColumns the name of the column, or the list of comma-separated columns
 * name composing the primary key of the referenced table.
 * @param onUpdate the action to do when a referenced primary key value is updated.
 * @param onDelete the action to do when a referenced record is deleted.
 */
void TableBuilder::addForeignKey(const char* constraintName, const char* fkColumns, const char* refTable, const char* refColumns, ReferenceOption onUpdate, ReferenceOption onDelete)
{
	ForeignKeyDesc desc;
	desc.constraintName = constraintName;
	desc.fkColumns = fkColumns;
	desc.refTable = refTable;
	desc.refColumns = refColumns;
	desc.onUpdate = onUpdate;
	desc.onDelete = onDelete;
	foreignKeys.push_back(desc);
}

/**
 * Adds a unicity constraint on a set of column to the table definition.
 * @param constraintName the name for the constraint.
 * @param columnsList name of the column whose values should be unique, or the
 * list of comma-separated columns name which should be unique tuples.
 */
void TableBuilder::addUnique(const char* constraintName, const char* columnsList)
{
	UniqueDesc desc;
	desc.constraintName = constraintName;
	desc.columnsList = columnsList;
	uniqueColumns.push_back(desc);
}

/**
 * Adds a expression checking constraint to the table definition.
 * @param constraintName the name for the constraint.
 * @param checkExpr an expression to be checked by the database each time a record
 * is updated or inserted. (example: "gender IN ('M','F')")
 */
void TableBuilder::addCheck(const char* constraintName, const char* checkExpr)
{
	CheckDesc desc;
	desc.constraintName = constraintName;
	desc.checkExpr = checkExpr;
	checks.push_back(desc);
}

/**
 * Generates the query using the initialization data provided by the previous other
 * method calls.
 */
std::string TableBuilder::buildQuery()
{
	std::ostringstream createQuery;
	beforeBuild(createQuery);

	createQuery << "CREATE TABLE " << dbName << '.' << tableName << "(\n";

	for(size_t i = 0; i < columns.size(); ++i)
	{
		const ColumnDesc& column = columns[i];
		if(i != 0)
		{
			createQuery << ",\n";
		}
		createQuery << '\t' << column.name << ' ' << getTypeName(column.type);
		if(column.length != UNDEFINED)
		{
			createQuery << '(' << column.length;
			if(column.decimals != UNDEFINED)
			{
				createQuery << ',' << column.decimals;
			}
			createQuery << ')';
		}
		if(column.sign == Unsigned)
		{
			createQuery << " UNSIGNED";
		}
		if(column.insertPolicy == AutoIncrement)
		{
			createQuery << " AUTO_INCREMENT";
		}
		else if(column.insertPolicy == CurrentDateTime)
		{
			createQuery << " DEFAULT " << getCurrentDateTimeSymbol();
		}
		if(column.nullPolicy == NotNull)
		{
			createQuery << " NOT NULL";
		}
	}
	if(primaryKeyColumns != 0)
	{
		createQuery << ",\n\t" <<
			"CONSTRAINT " << pkConstraintName << "\n\t\t" <<
			"PRIMARY KEY(" << primaryKeyColumns << ')';
	}
	for(size_t i = 0; i < foreignKeys.size(); ++i)
	{
		const ForeignKeyDesc& fkDesc = foreignKeys[i];
		createQuery << ",\n\t" <<
			"CONSTRAINT " << fkDesc.constraintName << "\n\t\t" << formatForeignKey(fkDesc);
		if(fkDesc.onDelete != NoAction)
		{
			createQuery << "\n\t\t"
				"ON DELETE " << getReferenceOptionName(fkDesc.onDelete);
		}
		if(fkDesc.onUpdate != NoAction)
		{
			createQuery << "\n\t\t"
				"ON UPDATE " << getReferenceOptionName(fkDesc.onUpdate);
		}
	}
	for(size_t i = 0; i < uniqueColumns.size(); ++i)
	{
		const UniqueDesc& ukDesc = uniqueColumns[i];
		createQuery << ",\n\t" <<
			"CONSTRAINT " << ukDesc.constraintName << "\n\t\t" <<
			"UNIQUE(" << ukDesc.columnsList << ')';
	}
	for(size_t i = 0; i < checks.size(); ++i)
	{
		const CheckDesc& ckDesc = checks[i];
		createQuery << ",\n\t" <<
			"CONSTRAINT " << ckDesc.constraintName << "\n\t\t" <<
			"CHECK(" << ckDesc.checkExpr << ')';
	}

	createQuery << std::endl << ')';
	afterBuild(createQuery);
	//std::cerr << "==DEBUG== TableBuilder::buildQuery\n" << createQuery.str() << std::endl;
	return createQuery.str();
}

/**
 * Return other queries to construct object that may be needed by the table.
 */
std::vector<std::string> TableBuilder::getSideQueries()
{
	return std::vector<std::string>();
}

/**
 * @fn const char* TableBuilder::getTypeName(ColumnType type) const
 * Returns the RDMS-specific type name used to represent this standard type concept.
 */

/**
 * Returns the RDMS-specific foreign key constraint.
 */
std::string TableBuilder::formatForeignKey(const ForeignKeyDesc & fkDesc)
{
	std::ostringstream fkConstraint;
	fkConstraint <<
			"FOREIGN KEY(" << fkDesc.fkColumns << ')' << "\n\t\t" <<
			"REFERENCES " << dbName << '.' << fkDesc.refTable <<
				'(' << fkDesc.refColumns << ')';
	return fkConstraint.str();
}

/**
 * @fn const char* TableBuilder::getReferenceOptionName(ReferenceOption option) const
 * Returns the RDBMS-specific name of the specified referential integrity option.
 */

/**
 * @fn const char* TableBuilder::getCurrentDateTimeSymbol() const
 * Return the name of the RDBMS-specific symbol that is interpreted as the current date and time.
 */

/**
 * This is a "callback" method called before building the query.
 * It can be re-implemented in RDBMS specific table builders to cope with its
 * particularities.
 */
void TableBuilder::beforeBuild(std::ostringstream& /*query*/)
{
}

/**
 * This is a "callback" method called after building the query.
 * It can be re-implemented in RDBMS specific table builders to cope with its
 * particularities.
 */
void TableBuilder::afterBuild(std::ostringstream& /* query */)
{
}


}
}
}
