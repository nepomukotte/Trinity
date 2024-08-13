#include "CCfg/Io/Db/OracleConfigDao.h"
#include "CCfg/Io/Db/ConfigBatch.h"
#include "CCfg/Io/Db/ConfigStructureVersionBatch.h"
#include "CCfg/Io/Db/ConfigVersionBatch.h"
#include "CCfg/Io/Db/ConfigClassBatch.h"
#include "CCfg/Io/Db/ConfigClassInfoBatch.h"
#include "CCfg/Io/Db/AbstractUpdateController.h"
#include "CCfg/Io/Db/Config.h"
#include "CCfg/Io/Db/Site.h"
#include "CCfg/Io/Db/SiteBatch.h"
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

struct OracleConfigDaoData
{
	static std::auto_ptr<OracleConfigDao> _dao;

	static OracleConfigDao* dao()
	{
		if(_dao.get() == 0)
		{
			_dao.reset(new OracleConfigDao);
			_dao->connect(TestParameters::instance().connectionString[CCfg::Io::Oracle]);
		}
		if(not _dao->isConnected())
		{
			throw std::runtime_error("DAO is not connected.");
		}
		return _dao.get();
	}
};

std::auto_ptr<OracleConfigDao> OracleConfigDaoData::_dao;

typedef test_group<OracleConfigDaoData> TestGroup;
typedef TestGroup::object TestObject;


///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

//#undef MYSQL_RDBMS
#define ORACLE_RDBMS
#include "ConfigDaoTest.hpp"

///////////////////////////////////////////////////////////////////////////////
//Test Group instanciation
///////////////////////////////////////////////////////////////////////////////

TestGroup tgOracleConfigDao("CCfg::Io::Db::OracleConfigDao");

}
