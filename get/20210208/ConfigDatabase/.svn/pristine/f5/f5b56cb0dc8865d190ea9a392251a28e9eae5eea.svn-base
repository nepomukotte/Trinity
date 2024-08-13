#ifndef CCFG_IO_DB_CONFIGCLASS_H_INCLUDED
#define CCFG_IO_DB_CONFIGCLASS_H_INCLUDED

#include "EntitySets.h"
#include "ConfigEntity.h"
#include <boost/intrusive_ptr.hpp>
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
 * This entity stands for an equivalence class of config attributes.
 *
 * @author fchateau
 */
struct ConfigClass : public ConfigEntity
{
	static const char TABLE_NAME[];
	static const unsigned NAME_MAX_LENGTH = 31;
	static const unsigned PATH_MAX_LENGTH = 255;
	enum ColumnsIndexes {ID = 0, FK_PARENT, FK_CREATION_STRUCT, NAME, PATH};

	UID fkParent;
	UID fkCreationStructure;
	std::string name;
	std::string path;

	ConfigClass();
	virtual ~ConfigClass();

	UID getId() const { return id; }
	const std::string& getPath() const;

	bool operator==(const ConfigClass& entity);
	bool operator!=(const ConfigClass& entity) { return not operator==(entity); }

private:
	ConfigClass(const ConfigClass& rvalue);
	ConfigClass& operator=(const ConfigClass& rvalue);
};

typedef boost::intrusive_ptr<ConfigClass> ConfigClassPtr;
typedef boost::multi_index_container<
	ConfigClassPtr,
	indexed_by<
		ordered_unique<tag<Id>, BOOST_MULTI_INDEX_CONST_MEM_FUN(ConfigClass, UID, getId)>,
		ordered_unique<tag<Path>, BOOST_MULTI_INDEX_CONST_MEM_FUN(ConfigClass, const std::string&, getPath)>
	>
> ConfigClassPtrSet;
typedef std::auto_ptr<ConfigClassPtrSet> ConfigClassPtrSetPtr;
typedef ConfigClassPtrSet::index<Id>::type ConfigClassIdIdx;
typedef ConfigClassPtrSet::index<Path>::type ConfigClassPathIdx;

}
}
}

#endif // CCFG_IO_DB_CONFIGCLASS_H_INCLUDED
