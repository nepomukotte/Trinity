#ifndef CCFG_IO_DB_CONFIGATTRIBUTE_H_INCLUDED
#define CCFG_IO_DB_CONFIGATTRIBUTE_H_INCLUDED

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
 * This entity stands for a scalar (variable) or an object (container) in a configuration.
 *
 * @author fchateau
 */
struct ConfigAttribute : public ConfigEntity
{
	static const char TABLE_NAME[];
	static const unsigned INDEX_MAX_LENGTH = 31;
	static const unsigned VALUE_MAX_LENGTH = 255;
	enum ColumnsIndexes {ID = 0, FK_PARENT, FK_VERSION, FK_CLASSINFO, INDEX, VALUE};

	UID fkParent;
	UID fkVersion;
	UID fkClassInfo;
	std::string index;
	std::string value;
	void* userData;

	ConfigAttribute();
	virtual ~ConfigAttribute();

	UID getId() const { return id; }

	bool operator==(const ConfigAttribute& entity);
	bool operator!=(const ConfigAttribute& entity) { return not operator==(entity); }

private:
	ConfigAttribute(const ConfigAttribute& rvalue);
	ConfigAttribute& operator=(const ConfigAttribute& rvalue);
};

typedef boost::intrusive_ptr<ConfigAttribute> ConfigAttributePtr;
typedef boost::multi_index_container<
	ConfigAttributePtr,
	indexed_by<
		ordered_unique<tag<Id>, BOOST_MULTI_INDEX_CONST_MEM_FUN(ConfigAttribute, UID, getId)>
	>
> ConfigAttributePtrSet;
typedef std::auto_ptr<ConfigAttributePtrSet> ConfigAttributePtrSetPtr;
typedef ConfigAttributePtrSet::index<Id>::type ConfigAttributeIdIdx;

}
}
}

#endif // CCFG_IO_DB_CONFIGATTRIBUTE_H_INCLUDED
