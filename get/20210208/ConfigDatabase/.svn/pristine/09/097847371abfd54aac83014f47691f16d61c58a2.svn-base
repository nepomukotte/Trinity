#include "ConfigClass.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigClass class methods body
////////////////////////////////////////////////////////////////////////////////

const char ConfigClass::TABLE_NAME[] = "CCfg_ConfigClass";
const unsigned ConfigClass::NAME_MAX_LENGTH;
const unsigned ConfigClass::PATH_MAX_LENGTH;

/**
 * Default constructor.
 */
ConfigClass::ConfigClass()
{
}

/**
 * Destructor.
 */
ConfigClass::~ConfigClass()
{
}

/**
 * Returns the path of the configuration class.
 */
const std::string& ConfigClass::getPath() const
{
	return path;
}

/**
 * Equality operator.
 */
bool ConfigClass::operator==(const ConfigClass& entity)
{
	return id == entity.id and
		fkParent == entity.fkParent and
		fkCreationStructure == entity.fkCreationStructure and
		name == entity.name and
		path == entity.path;
}

/**
 * @fn bool ConfigClass::operator!=(const ConfigClass& entity)
 * Unequality operator.
 */

}
}
}
