#ifndef CCFG_IO_DB_OBJRELMAPPING_H_INCLUDED
#define CCFG_IO_DB_OBJRELMAPPING_H_INCLUDED

#include "Config.h"
#include "ConfigAttribute.h"
#include "ConfigClass.h"
#include "ConfigClassInfo.h"
#include "ConfigDatabase.h"
#include "ConfigStructureVersion.h"
#include "ConfigVersion.h"
#include <soci.h>


namespace SOCI
{
using namespace CCfg::Io::Db;

template<>
struct TypeConversion<ConfigDatabasePtr>
{
	typedef Values base_type;
	static ConfigDatabasePtr from(const Values& v)
	{
		ConfigDatabasePtr db(new ConfigDatabase);
		db->name = v.get<std::string>(0);
		return db;
	}

	static Values to(ConfigDatabasePtr& db)
	{
		Values v;
		v.set("schema_name", db->name);
		return v;
	}
};

}

#endif // CCFG_IO_DB_OBJRELMAPPING_H_INCLUDED
