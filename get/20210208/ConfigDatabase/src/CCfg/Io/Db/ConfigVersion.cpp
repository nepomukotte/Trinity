#include "ExtractionProcess.h"
#include "ConfigVersion.h"
#include "ConfigDao.h"
#include "../../../Utils/Date/StdTmHelpers.h"
#include <cassert>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigVersion class methods body
////////////////////////////////////////////////////////////////////////////////

const char ConfigVersion::TABLE_NAME[] = "CCfg_ConfigVersion";
const unsigned ConfigVersion::NAME_MAX_LENGTH;

/**
 * Default constructor.
 */
ConfigVersion::ConfigVersion()
{
}

/**
 * Destructor.
 */
ConfigVersion::~ConfigVersion()
{
}

/**
 * Returns the date of creation of this version.
 */
const std::tm& ConfigVersion::getCreationDate() const
{
	return creationDate;
}

/**
 * Return the date of creation of this version as a time_t (number of seconds
 * since the Epoch).
 */
time_t ConfigVersion::getCreationDateAsTimeT()
{
	return mktime(&creationDate);
}

/**
 * Equality operator.
 */
bool ConfigVersion::operator==(const ConfigVersion& entity)
{
	return id == entity.id and
		fkStructure == entity.fkStructure and
		creationDate == entity.creationDate and
		name == entity.name;
}

/**
 * @fn bool ConfigVersion::operator!=(const ConfigVersion& entity)
 * Unequality operator.
 */

/**
 * Extracts a configuration from the database.
 * @return Returns the created configuration tree.
 */
std::auto_ptr<Attribute> ConfigVersion::extractConfig()
{
	ExtractionProcess extractor(*this);
	return extractor.extract();
}

/**
 * Retrieves the structure version associated with this version, from the database.
 */
ConfigStructureVersionPtr ConfigVersion::getStructure()
{
	assert(dao != 0);
	return dao->getStructure(fkStructure);
}

/**
 * Retrieves the list of attributes associated with this version, from the database.
 */
ConfigAttributePtrSetPtr ConfigVersion::getAttributes()
{
	assert(dao != 0);
	return dao->getAttributesOfVersion(*this);
}

}
}
}
