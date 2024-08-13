#include "CCfg/Io/Db/PostgreSqlConfigDao.h"
#include "CCfg/Io/Db/ConfigBatch.h"
#include "CCfg/Io/Db/ConfigStructureVersionBatch.h"
#include "CCfg/Io/Db/ConfigVersionBatch.h"
#include "CCfg/Io/Db/ConfigClassBatch.h"
#include "CCfg/Io/Db/ConfigClassInfoBatch.h"
#include "CCfg/Io/Db/Config.h"
#include "CCfg/Io/Db/Site.h"
#include "CCfg/Io/Db/SiteBatch.h"
#include "CCfg/Io/Db/AbstractUpdateController.h"
#include "CCfg/Io/Document.h"
#include "CCfg/StructureNode.h"
#include "CCfg/Object.h"
#include "CCfg/Io/DatabaseStorage.h"
#include "Utils/Date/StdTmHelpers.h"
#include "../../../TestParameters.h"
#include "Sleep.h"
#include "tut.h"
#include <memory>

using namespace CCfg::Io::Db;

namespace tut
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct PostgreSqlConfigDaoData
{
	static std::auto_ptr<PostgreSqlConfigDao> _dao;

	static PostgreSqlConfigDao* dao()
	{
		if(_dao.get() == 0)
		{
			_dao.reset(new PostgreSqlConfigDao);
			_dao->connect(TestParameters::instance().connectionString[CCfg::Io::PostgreSql]);
		}
		if(not _dao->isConnected())
		{
			throw std::runtime_error("DAO is not connected.");
		}
		return _dao.get();
	}
};

std::auto_ptr<PostgreSqlConfigDao> PostgreSqlConfigDaoData::_dao;

typedef test_group<PostgreSqlConfigDaoData> TestGroup;
typedef TestGroup::object TestObject;


///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

//#undef ORACLE_RDBMS
#define POSTGRESQL_RDBMS
#include "ConfigDaoTest.hpp"


///////////////////////////////////////////////////////////////////////////////
//Test Group instantiation
///////////////////////////////////////////////////////////////////////////////

TestGroup tgPostgreSqlConfigDao("CCfg::Io::Db::PostgreSqlConfigDao");

}
