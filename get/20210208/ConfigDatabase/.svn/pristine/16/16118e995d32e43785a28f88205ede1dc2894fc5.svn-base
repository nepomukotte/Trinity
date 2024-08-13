#ifndef CCFG_IO_DB_CONFIGVERSION_H_INCLUDED
#define CCFG_IO_DB_CONFIGVERSION_H_INCLUDED

#include "EntitySets.h"
#include "ConfigAttribute.h"
#include "ConfigEntity.h"
#include "CCfg/Attribute.h"
#include <boost/intrusive_ptr.hpp>
#include <string>
#include <ctime>

namespace CCfg
{
namespace Io
{
namespace Db
{
// forward declarations
struct ConfigStructureVersion;
typedef boost::intrusive_ptr<ConfigStructureVersion> ConfigStructureVersionPtr;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This entity stands for a version of the values of a configuration.
 *
 * @author fchateau
 */
struct ConfigVersion : public ConfigEntity
{
	static const char TABLE_NAME[];
	static const unsigned NAME_MAX_LENGTH = 255;
	enum ColumnsIndexes {ID = 0, CREATION_DATE, NAME, FK_STRUCTURE};

	UID fkStructure;
	std::tm creationDate;
	std::string name;

	ConfigVersion();
	virtual ~ConfigVersion();

	UID getId() const { return id; }
	const std::tm& getCreationDate() const;
	time_t getCreationDateAsTimeT();

	bool operator==(const ConfigVersion& entity);
	bool operator!=(const ConfigVersion& entity) { return not operator==(entity); }

	std::auto_ptr<Attribute> extractConfig();
	ConfigStructureVersionPtr getStructure();
	ConfigAttributePtrSetPtr getAttributes();

private:
	ConfigVersion(const ConfigVersion& rvalue);
	ConfigVersion& operator=(const ConfigVersion& rvalue);
};

typedef boost::intrusive_ptr<ConfigVersion> ConfigVersionPtr;
typedef boost::multi_index_container<
	ConfigVersionPtr,
	indexed_by<
		ordered_unique<tag<Id>, BOOST_MULTI_INDEX_CONST_MEM_FUN(ConfigVersion, UID, getId)>,
		ordered_unique<tag<Date>, BOOST_MULTI_INDEX_MEM_FUN(ConfigVersion, time_t, getCreationDateAsTimeT)>
	>
> ConfigVersionPtrSet;
typedef std::auto_ptr<ConfigVersionPtrSet> ConfigVersionPtrSetPtr;
typedef ConfigVersionPtrSet::index<Id>::type ConfigVersionIdIdx;
typedef ConfigVersionPtrSet::index<Date>::type ConfigVersionDateIdx;

}
}
}

#endif // CCFG_IO_DB_CONFIGVERSION_H_INCLUDED
