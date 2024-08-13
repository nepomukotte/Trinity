#include "CCfg/Io/Db/Sqlite3ConfigDao.h"
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

struct Sqlite3ConfigDaoData
{
	static std::auto_ptr<Sqlite3ConfigDao> _dao;

	static Sqlite3ConfigDao* dao()
	{
		if(_dao.get() == 0)
		{
			_dao.reset(new Sqlite3ConfigDao);
			_dao->connect(TestParameters::instance().connectionString[CCfg::Io::Sqlite3]);
		}
		if(not _dao->isConnected())
		{
			throw std::runtime_error("DAO is not connected.");
		}
		return _dao.get();
	}
};

std::auto_ptr<Sqlite3ConfigDao> Sqlite3ConfigDaoData::_dao;

typedef test_group<Sqlite3ConfigDaoData> TestGroup;
typedef TestGroup::object TestObject;


///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

#define SQLITE3_RDBMS
#include "ConfigDaoTest.hpp"


///////////////////////////////////////////////////////////////////////////////
//Test Group instanciation
///////////////////////////////////////////////////////////////////////////////

TestGroup tgSqlite3ConfigDao("CCfg::Io::Db::Sqlite3ConfigDao");

}
