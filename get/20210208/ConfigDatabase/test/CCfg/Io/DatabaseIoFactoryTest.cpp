#include "Db/Sleep.h"
#include "CCfg/Io/Db/Config.h"
#include "CCfg/Io/Db/ConfigVersion.h"
#include "CCfg/Io/Db/ConfigDao.h"
#include "CCfg/Io/FormatFactory.h"
#include "CCfg/Io/StorageFactory.h"
#include "CCfg/Io/DatabaseFormatFactory.h"
#include "CCfg/Io/DatabaseStorageFactory.h"
#include "CCfg/Io/UnhandledUrl.h"
#include "../../TestParameters.h"
#include "CCfg/Io/Document.h"
#include "CCfg/Io/Url.h"
#include "CCfg/Object.h"
#include "tut.h"

namespace tut
{
using namespace CCfg::Io;
using namespace CCfg::Io::Db;

////////////////////////////////////////////////////////////////////////////////
// a few classes for tests
////////////////////////////////////////////////////////////////////////////////

struct ValidLoginCallback : AbstractLoginCallback
{
	bool requestLogin(DbmsType dbms,
			const std::string& host, const std::string& port,
			const std::string& sid, std::string& user, std::string& pwd)
	{
		ensure_equals("host", host, TestParameters::instance().host[dbms]);
		ensure_equals("port", port, TestParameters::instance().port[dbms]);
		ensure_equals("user", user, TestParameters::instance().user[dbms]);
		ensure_equals("sid", sid, TestParameters::instance().sid[dbms]);
		pwd = TestParameters::instance().passwd[dbms];
		return true;
	}
};

/**
 * This class is used to registers the database format factory exactly once.
 */
class FormatRegisterer
{
public:
	~FormatRegisterer() {}

	static FormatRegisterer& instance()
	{
		static FormatRegisterer inst;
		return inst;
	}

private:
	FormatRegisterer()
	{
		DatabaseFormatFactory* frFactory = new DatabaseFormatFactory;
		std::auto_ptr<AbstractFormatFactory> factory1(frFactory);
		frFactory->setUpdateController(UpdateControllerPtr(new Db::DefaultUpdateController));
		FormatFactory::instance().registerFormat(factory1);

		DatabaseStorageFactory* stFactory = new DatabaseStorageFactory;
		std::auto_ptr<AbstractStorageFactory> factory2(stFactory);
		stFactory->setLoginCallback(LoginCallbackPtr(new ValidLoginCallback));
		StorageFactory::instance().registerStorage(factory2);
	}
};


////////////////////////////////////////////////////////////////////////////////
// test fixture definition
////////////////////////////////////////////////////////////////////////////////

struct DbIoFactoryTest
{
	DbIoFactoryTest()
	{
		FormatRegisterer::instance();
	}
	~DbIoFactoryTest()
	{
	}

	void testIoFactory(DbmsType dbms);
};

typedef test_group<DbIoFactoryTest> testgroup;
typedef testgroup::object testobject;


////////////////////////////////////////////////////////////////////////////////
// test methods body
////////////////////////////////////////////////////////////////////////////////

void registerFactory()
{

}

/**
 * Test configuration insertion and extraction for the specified database.
 */
void DbIoFactoryTest::testIoFactory(DbmsType dbms)
{
	// build the url
	Url url;
	DatabaseStorage::fillUrl(url, dbms,
			TestParameters::instance().user[dbms],
			std::string(), // clear password: it will be set by the login callback
			TestParameters::instance().host[dbms],
			TestParameters::instance().port[dbms],
			TestParameters::instance().sid[dbms],
			TestParameters::instance().schema[dbms],
			"config1");

	// create storage and format object
	std::auto_ptr<Format> format = FormatFactory::instance().createFormat(url);
	std::auto_ptr<Storage> storage = StorageFactory::instance().createStorage(url);
	DatabaseStorage& dbStorage = dynamic_cast<DatabaseStorage&>(*storage);
	// open the storage and reinitialize it
	dbStorage.open();
	Db::ConfigDao& dao = dbStorage.getDao();
	dao.replaceAllTables();

	// check that the login callback set the password
	ensure_equals("1", dbStorage.getPassword(), TestParameters::instance().passwd[dbms]);

	// load a configuration from a file
	Document doc1;
	CCfg::Object* cfg = &doc1.loadFromFile("testconfig1.xml");
	// insert the configuration into the database
	format->format(cfg, *storage, NoOpIncludeHandler::getInstance());
	// extraction the configuration
	std::auto_ptr<CCfg::Object> dbCfg = format->parse(*storage, NoOpIncludeHandler::getInstance());
	// compares it with the original
	ensure("2", cfg->equals(*dbCfg, CCfg::EQUIVALENT));

	// load a second configuration from a file
	Document doc2;
	CCfg::Object* cfg2 = &doc2.loadFromFile("testconfig2.xml");
	// insert the configuration as an update
	sleep(1); // avoid timestamp doublets
	format->format(cfg2, *storage, NoOpIncludeHandler::getInstance());

	// retrieve the versions records of the 2 configuration we have just inserted.
	ConfigPtr config1 = dao.searchConfig("config1");
	ensure("3", config1.get() != 0);
	ConfigVersionPtrSetPtr versions = config1->getVersions();
	ensure("4", versions.get() != 0);
	ensure_equals("5", versions->size(), 2u);

	// retrieve the first configuration
	ConfigVersionIdIdx::iterator curVersion = versions->get<Id>().begin();
	dbStorage.setCurrentConfigDate((**curVersion).creationDate);
	dbCfg = format->parse(*storage, NoOpIncludeHandler::getInstance());
	ensure("6", cfg->equals(*dbCfg, CCfg::EQUIVALENT));
	// retrieve the second configuration and check it
	++curVersion;
	dbStorage.setCurrentConfigDate((**curVersion).creationDate);
	dbCfg = format->parse(*storage, NoOpIncludeHandler::getInstance());
	ensure("7", cfg2->equals(*dbCfg, CCfg::EQUIVALENT));
}

/**
 * Test URL retrieval.
 */
template<> template<>
void testobject::test<1>()
{
	std::tm date;
	date.tm_year = 100;
	date.tm_mon = 0;
	date.tm_mday = 0;
	date.tm_hour = 12;
	date.tm_min = 5;
	date.tm_sec = 10;
	Url url;
	DatabaseStorage::fillUrl(url, Oracle,
			TestParameters::instance().user[Oracle],
			std::string(),
			TestParameters::instance().host[Oracle],
			TestParameters::instance().port[Oracle],
			TestParameters::instance().sid[Oracle],
			TestParameters::instance().schema[Oracle],
			"config1", &date);
	std::auto_ptr<Storage> storage = StorageFactory::instance().createStorage(url);
	const std::string& urlStr = storage->getUrl();
	Url url2(urlStr);
	ensure(url == url2);
}

/**
 * Test storage construction.
 */
template<> template<>
void testobject::test<2>()
{
	Url url;
	DatabaseStorage::fillUrl(url, Oracle,
				TestParameters::instance().user[Oracle],
				TestParameters::instance().passwd[Oracle],
				TestParameters::instance().host[Oracle],
				TestParameters::instance().port[Oracle],
				TestParameters::instance().sid[Oracle],
				TestParameters::instance().schema[Oracle],
				"config1");
	std::auto_ptr<Storage> storage = StorageFactory::instance().createStorage(url);
	DatabaseStorage& dbStorage = dynamic_cast<DatabaseStorage&>(*storage);

	ensure_equals(dbStorage.getDbmsType(), Oracle);
	ensure_equals(dbStorage.getHost(), TestParameters::instance().host[Oracle]);
	ensure_equals(dbStorage.getPort(), TestParameters::instance().port[Oracle]);
	ensure_equals(dbStorage.getUser(), TestParameters::instance().user[Oracle]);
	ensure_equals(dbStorage.getPassword(), TestParameters::instance().passwd[Oracle]);
	ensure_equals(dbStorage.getSID(), TestParameters::instance().sid[Oracle]);
	ensure_equals(dbStorage.getSchema(), TestParameters::instance().schema[Oracle]);
	ensure_equals(dbStorage.getCurrentConfigName(), "config1");
}

#ifdef CCFG_USE_ORACLE
/**
 * Test configuration insertion and extraction through Oracle databases.
 */
template<> template<>
void testobject::test<3>()
{
	testIoFactory(Oracle);
}
#endif

/**
 * Test configuration insertion and extraction through MySql databases.
 */
template<> template<>
void testobject::test<4>()
{
	testIoFactory(MySql);
}

#ifdef CCFG_USE_POSTGRESQL
/**
 * Test configuration insertion and extraction through PostgreSql databases.
 */
template<> template<>
void testobject::test<5>()
{
	testIoFactory(PostgreSql);
}
#endif

/**
 * Test what happens when we use an invalid database URL.
 */
template<> template<>
void testobject::test<6>()
{
	try
	{
		Url url("foo:///home/user/foo.xml"); // use an invalid url for database IO
		std::auto_ptr<Format> format = FormatFactory::instance().createFormat(url);
		fail("The Url should produce an error.");
	}
	catch(const UnhandledUrl& e)
	{
		ensure_equals(e.getUrl(), "foo:///home/user/foo.xml");
		ensure_equals(e.getProduct(), "Format");
	}
}

/**
 * Test what happens when trying to connect to an unexisting database.
 */
template<> template<>
void testobject::test<7>()
{
	try
	{
		Url url("mysql://foo:42@bar:1521/foo/config1");
		std::auto_ptr<Format> format = FormatFactory::instance().createFormat(url);
		std::auto_ptr<Storage> storage = StorageFactory::instance().createStorage(url);
		DatabaseStorage& dbStorage = dynamic_cast<DatabaseStorage&>(*storage);
		dbStorage.open();
		fail("The connection should fail");
	}
	catch(const DatabaseStorage::ConnectionError& e)
	{
	}
}

////////////////////////////////////////////////////////////////////////////////
//Test Group instantiation
////////////////////////////////////////////////////////////////////////////////

testgroup tgDatabaseIoFactory("CCfg::Io::DatabaseIoFactory");

}
