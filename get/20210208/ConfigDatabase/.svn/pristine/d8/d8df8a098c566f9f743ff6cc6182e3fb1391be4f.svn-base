#ifndef CCFG_IO_DB_CONFIGCLASSINFO_H_INCLUDED
#define CCFG_IO_DB_CONFIGCLASSINFO_H_INCLUDED

#include "EntitySets.h"
#include "ConfigEntity.h"
#include "CCfg/TypeInfo.h"
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
 * This entity stands for the informations concerning equivalence classes that
 * may vary when the structure changes.
 *
 * @author fchateau
 */
struct ConfigClassInfo : public ConfigEntity
{
	static const char TABLE_NAME[];
	static const unsigned UNIT_MAX_LENGTH = 31;
	static const unsigned RANGE_MAX_LENGTH = 255;
	enum ColumnsIndexes {ID = 0, TYPE, UNIT, RANGE, FK_CLASS, FK_STRUCTURE};

	UID fkClass;
	UID fkStructure;
	AttributeTypeId type;
	std::string unit;
	std::string range;

	ConfigClassInfo();
	virtual ~ConfigClassInfo();

	UID getId() const { return id; }
	UID getFkClass() const;

	bool operator==(const ConfigClassInfo& entity);
	bool operator!=(const ConfigClassInfo& entity) { return not operator==(entity); }

private:
	ConfigClassInfo(const ConfigClassInfo& rvalue);
	ConfigClassInfo& operator=(const ConfigClassInfo& rvalue);
};

typedef boost::intrusive_ptr<ConfigClassInfo> ConfigClassInfoPtr;
typedef boost::multi_index_container<
	ConfigClassInfoPtr,
	indexed_by<
		ordered_unique<tag<Id>, BOOST_MULTI_INDEX_CONST_MEM_FUN(ConfigClassInfo, UID, getId)>,
		ordered_non_unique<tag<Cfgclass>, BOOST_MULTI_INDEX_CONST_MEM_FUN(ConfigClassInfo, UID, getFkClass)>
	>
> ConfigClassInfoPtrSet;
typedef std::auto_ptr<ConfigClassInfoPtrSet> ConfigClassInfoPtrSetPtr;
typedef ConfigClassInfoPtrSet::index<Id>::type ConfigClassInfoIdIdx;
typedef ConfigClassInfoPtrSet::index<Cfgclass>::type ConfigClassInfoClassIdx;

}
}
}

#endif // CCFG_IO_DB_CONFIGCLASSINFO_H_INCLUDED
