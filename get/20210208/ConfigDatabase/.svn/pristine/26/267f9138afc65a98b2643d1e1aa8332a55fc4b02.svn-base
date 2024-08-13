
/**
 * Test table creation.
 */
template<> template<>
void TestObject::test<1>()
{
	dao()->replaceAllTables();
}

/**
 * Test database enumeration.
 */
template<> template<>
void TestObject::test<2>()
{
	ConfigDatabasePtrSetPtr databases = dao()->getDatabases();
	ConfigDatabaseNameIdx& nameIdx = databases->get<Name>();
	#ifdef MYSQL_RDBMS
		ensure(nameIdx.find("information_schema") != nameIdx.end());
	#endif
	#ifdef ORACLE_RDBMS
		ensure(nameIdx.find("SYSTEM") != nameIdx.end());
		ensure(nameIdx.find("SYS") != nameIdx.end());
	#endif
	#ifdef POSTGRESQL_RDBMS
		ensure(true);
	#endif
	#ifdef SQLITE3_RDBMS
		ensure(nameIdx.find("main") != nameIdx.end());
	#endif
}

/**
 * Test Config insertion and selection.
 */
template<> template<>
void TestObject::test<3>()
{
	// create a batch of 3 Config objects.
	ConfigBatch batch;
	Config config;

	UID uid = dao()->allocIdRangeForConfigs(3u);

	config.id = ++uid;
	config.name = "config1";
	batch.append(config);

	config.id = ++uid;
	config.name = "config2";
	batch.append(config);

	config.id = ++uid;
	config.name = "config3";
	batch.append(config);

	// insert the 3 previous configurations.
	dao()->insert(batch);

	// retrieve the inserted configurations.
	ConfigPtrSetPtr configs = dao()->getConfigs();

	// check if we retrieved the right number of configurations.
	ensure_equals(configs->size(), 3u);

	// iterate on retrieved configurations to compare them with those in the batch.
	ConfigIdIdx& idIdx = configs->get<Id>();
	UID i = 0;
	for(ConfigIdIdx::const_iterator it = idIdx.begin(); it != idIdx.end(); ++it, ++i)
	{
		batch.assignEntity(i, config);
		ensure(**it == config);
	}

	// retrieve the next id
	uid = dao()->allocIdRangeForConfigs(0);

	// check if it has been incremented correctly
	ensure_equals(uid, 3u);
}

/**
 * Test ConfigStructureVersion insertion and selection.
 */
template<> template<>
void TestObject::test<4>()
{
	// create a batch of 5 ConfigStructureVersion objects.
	ConfigStructureVersionBatch batch;
	ConfigStructureVersion sversion;

	UID uid = dao()->allocIdRangeForStructureVersions(5u);

	// the first 3 versions are related to config #1
	sversion.id = ++uid;
	sversion.creationDate.tm_sec = 0;
	sversion.creationDate.tm_min = 0;
	sversion.creationDate.tm_hour = 0;
	sversion.creationDate.tm_mday = 1;
	sversion.creationDate.tm_mon = 0;
	sversion.creationDate.tm_year = 100;
	sversion.fkConfig = 1;
	batch.append(sversion);

	sversion.id = ++uid;
	sversion.creationDate.tm_mon = 1;
	batch.append(sversion);

	sversion.id = ++uid;
	sversion.creationDate.tm_mon = 2;
	batch.append(sversion);

	// the last 2 versions are related to config #2
	sversion.id = ++uid;
	sversion.creationDate.tm_mon = 0;
	sversion.creationDate.tm_year = 105;
	sversion.fkConfig = 2;
	batch.append(sversion);

	sversion.id = ++uid;
	sversion.creationDate.tm_mday = 2;
	batch.append(sversion);

	// insert the 5 ConfigStructureVersion.
	dao()->insert(batch);

	// retrieve the inserted ConfigStructureVersion.
	ConfigPtr config = dao()->getConfig(1);
	ConfigStructureVersionPtrSetPtr sversionSet = config->getStructureVersions();
	ensure_equals(sversionSet->size(), 3u);

	// iterate on retrieved ConfigStructureVersion to compare them with those in the batch.
	ConfigStructureVersionIdIdx& idIdx = sversionSet->get<Id>();
	UID i = 0;
	for(ConfigStructureVersionIdIdx::const_iterator it = idIdx.begin(); it != idIdx.end(); ++it, ++i)
	{
		batch.assignEntity(i, sversion);
		ensure(**it == sversion);
	}

	// retrieve the next id
	uid = dao()->allocIdRangeForStructureVersions(0);

	// check if it has been incremented correctly
	ensure_equals(uid, 5u);
}

/**
 *
 */
template<> template<>
void TestObject::test<5>()
{
	// create a batch of 3 ConfigVersion objects.
	ConfigVersionBatch batch;
	ConfigVersion version;

	UID uid = dao()->allocIdRangeForVersions(3u);

	// create 1 ConfigVersion per ConfigStructureVersion for Config #1.
	version.id = ++uid;
	version.creationDate.tm_sec = 0;
	version.creationDate.tm_min = 0;
	version.creationDate.tm_hour = 0;
	version.creationDate.tm_mday = 1;
	version.creationDate.tm_mon = 0;
	version.creationDate.tm_year = 100;
	version.fkStructure = 1;
	version.name = "version1";
	batch.append(version);

	version.id = ++uid;
	version.creationDate.tm_mon = 1;
	version.fkStructure = 2;
	version.name = "version2";
	batch.append(version);

	version.id = ++uid;
	version.creationDate.tm_mon = 2;
	version.fkStructure = 3;
	version.name = "version3";
	batch.append(version);

	// insert the 3 ConfigVersion.
	dao()->insert(batch);

	// retrieve the inserted ConfigVersion.
	ConfigPtr config = dao()->getConfig(1);
	ConfigVersionPtrSetPtr versionSet = config->getVersions();
	ensure_equals(versionSet->size(), 3u);

	// iterate on retrieved ConfigVersion to compare them with those in the batch.
	ConfigVersionIdIdx& idIdx = versionSet->get<Id>();
	UID i = 0;
	for(ConfigVersionIdIdx::const_iterator it = idIdx.begin(); it != idIdx.end(); ++it, ++i)
	{
		batch.assignEntity(i, version);
		ensure(**it == version);
	}

	// retrieve the next id
	uid = dao()->allocIdRangeForVersions(0);

	// check if it has been incremented correctly
	ensure_equals(uid, 3u);
}

/**
 *
 */
template<> template<>
void TestObject::test<6>()
{
	// create a batch of 3 ConfigClass objects.
	ConfigClassBatch batch;
	ConfigClass clazz;

	UID uid = dao()->allocIdRangeForClasses(3u);

	// create 3 classes for the first version
	clazz.id = ++uid;
	clazz.fkParent = NULL_UID;
	clazz.fkCreationStructure = 1;
	clazz.name = "A";
	clazz.path = "A";
	batch.append(clazz);

	clazz.id = ++uid;
	clazz.fkParent = 1;
	clazz.fkCreationStructure = 1;
	clazz.name = "x";
	clazz.path = "A.x";
	batch.append(clazz);

	clazz.id = ++uid;
	clazz.fkParent = 1;
	clazz.fkCreationStructure = 1;
	clazz.name = "y";
	clazz.path = "A.y";
	batch.append(clazz);

	// insert the 3 ConfigClass.
	dao()->insert(batch);

	// create a class info for each class
	// this is necessary because getClassesOfStructure returns only classes for
	// which a class info exists for the structure.
	ConfigClassInfoBatch batch2;
	ConfigClassInfo info;

	uid = dao()->allocIdRangeForClassInfos(3u);

	info.id = ++uid;
	info.fkClass = 1;
	info.fkStructure = 1;
	info.type = CCfg::OBJECT;
	info.unit = "";
	info.range = "";
	batch2.append(info);

	info.id = ++uid;
	info.fkClass = 2;
	info.fkStructure = 1;
	info.type = CCfg::INT_SCALAR;
	info.unit = "mm";
	info.range = "[0,1000]";
	batch2.append(info);

	info.id = ++uid;
	info.fkClass = 3;
	info.fkStructure = 1;
	info.type = CCfg::INT_SCALAR;
	info.unit = "mm";
	info.range = "[0,1000]";
	batch2.append(info);

	// insert the 3 ConfigClassInfo.
	dao()->insert(batch2);

	ConfigStructureVersionPtr structure = dao()->getStructure(1);

	// retrieve the inserted ConfigClass.
	ConfigClassPtrSetPtr classSet = structure->getClasses();
	ensure_equals(classSet->size(), 3u);

	// retrieve the inserted ConfigClassInfo.
	ConfigClassInfoPtrSetPtr infoSet = structure->getClassInfos();
	ensure_equals(infoSet->size(), 3u);

	// iterate on retrieved ConfigClass to compare them with those in the batch.
	ConfigClassIdIdx& idIdx = classSet->get<Id>();
	UID i = 0;
	for(ConfigClassIdIdx::const_iterator it = idIdx.begin(); it != idIdx.end(); ++it, ++i)
	{
		batch.assignEntity(i, clazz);
		ensure(**it == clazz);
	}

	// iterate on retrieved ConfigClassInfo to compare them with those in the batch.
	ConfigClassInfoIdIdx& idIdx2 = infoSet->get<Id>();
	UID j = 0;
	for(ConfigClassInfoIdIdx::const_iterator it = idIdx2.begin(); it != idIdx2.end(); ++it, ++j)
	{
		batch2.assignEntity(j, info);
		ensure(**it == info);
	}

	// retrieve the next ConfigClass id
	uid = dao()->allocIdRangeForClasses(0);

	// check if it has been incremented correctly
	ensure_equals(uid, 3u);

	// retrieve the next ConfigClassInfo id
	uid = dao()->allocIdRangeForClassInfos(0);

	// check if it has been incremented correctly
	ensure_equals(uid, 3u);
}

/**
 * Test method getLastStructureVersion
 */
template<> template<>
void TestObject::test<8>()
{
	// retrieve the last structure version of config1 and check its properties
	ConfigPtr config = dao()->searchConfig("config1");
	ConfigStructureVersionPtr structVersion = config->getLastStructureVersion();
	ensure(structVersion.get() != 0);
	ensure_equals(structVersion->id, 3u);
	ensure_equals(structVersion->fkConfig, 1u);
	ensure_equals(structVersion->creationDate.tm_mon, 2);

	// check that the last structure version of config3 does not exist
	config = dao()->searchConfig("config3");
	structVersion = config->getLastStructureVersion();
	ensure(structVersion.get() == 0);
}

/**
 * Test method getLastVersion
 */
template<> template<>
void TestObject::test<9>()
{
	// retrieve the last version of config1 and check its properties
	ConfigPtr config = dao()->searchConfig("config1");
	ConfigVersionPtr version = config->getLastVersion();
	ensure(version.get() != 0);
	ensure_equals(version->id, 3u);
	ensure_equals(version->fkStructure, 3u);
	ensure_equals(version->creationDate.tm_mon, 2);
	ensure_equals(version->name, "version3");

	// check that the last version of config3 does not exist
	config = dao()->searchConfig("config3");
	version = config->getLastVersion();
	ensure(version.get() == 0);
}

/**
 *
 */
template<> template<>
void TestObject::test<10>()
{
	CCfg::Io::Document doc;
	CCfg::Attribute& cfg = doc.loadFromFile("testconfig1.xml");
	struct TestUpdateController : AbstractUpdateController
	{
		TestUpdateController() {}
		~TestUpdateController() {}
		void onClassAdded(const CCfg::StructureNode& /* newClass */)
		{
			fail();
		}
		void onClassRemoved(const CCfg::StructureNode& /* oldClass */)
		{
			fail();
		}
		void onUnitChanged(const CCfg::StructureNode& /* newClass */, const CCfg::StructureNode& /* oldClass */)
		{
			fail();
		}
		void onRangeChanged(const CCfg::StructureNode& /* newClass */, const CCfg::StructureNode& /* oldClass */)
		{
			fail();
		}
		void onTypeChanged(const CCfg::StructureNode& /* newClass */, const CCfg::StructureNode& /* oldClass */)
		{
			fail();
		}
		bool doConfirmStructureChange()
		{
			fail();
			return false;
		}
	} control;
	ConfigPtr config = dao()->searchConfig("config3");
	ConfigVersionPtr version = config->insertNewVersion(&cfg, control);

	std::auto_ptr<CCfg::Attribute> cfg2 = version->extractConfig();
	ensure(cfg2->equals(cfg, CCfg::EQUIVALENT));
}

/**
 *
 */
template<> template<>
void TestObject::test<11>()
{
	CCfg::Io::Document doc;
	CCfg::Attribute& cfg = doc.loadFromFile("testconfig2.xml");
	struct TestUpdateController : AbstractUpdateController
	{
		TestUpdateController() {}
		~TestUpdateController() {}
		void onClassAdded(const CCfg::StructureNode& newClass)
		{
			ensure_equals(newClass.path, "cfg.Computer.HardDisk.RotationSpeed");
		}
		void onClassRemoved(const CCfg::StructureNode& oldClass)
		{
			ensure_equals(oldClass.path, "cfg.Computer.CPU.HyperThreading");
		}
		void onUnitChanged(const CCfg::StructureNode& newClass, const CCfg::StructureNode& oldClass)
		{
			ensure_equals(oldClass.path, "cfg.Computer.CPU.Frequency");
			ensure_equals(oldClass.unit, "MHz");
			ensure_equals(newClass.unit, "GHz");
		}
		void onRangeChanged(const CCfg::StructureNode& newClass, const CCfg::StructureNode& oldClass)
		{
			ensure_equals(oldClass.path, "cfg.Computer.HardDisk.Capacity");
			ensure_equals(newClass.range, "[1, 2000]");
		}
		void onTypeChanged(const CCfg::StructureNode& newClass, const CCfg::StructureNode& oldClass)
		{
			ensure_equals(oldClass.path, "cfg.Computer.CPU.Frequency");
			ensure_equals(oldClass.type, CCfg::INT_SCALAR);
			ensure_equals(newClass.type, CCfg::REAL_SCALAR);
		}
		bool doConfirmStructureChange()
		{
			return true;
		}
	} control;
	sleep(1); // wait 1 sec to avoid using the same date than the previous version.

	ConfigPtr config = dao()->searchConfig("config3");
	ConfigVersionPtr version = config->insertNewVersion(&cfg, control);

	std::auto_ptr<CCfg::Attribute> cfg2 = version->extractConfig();
	ensure(cfg2->equals(cfg, CCfg::EQUIVALENT));
}

/**
 *
 */
template<> template<>
void TestObject::test<12>()
{
	ConfigPtr config = dao()->searchConfig("config999");
	ensure(config.get() == 0);
}

/**
 * Test site insertion.
 */
template<> template<>
void TestObject::test<13>()
{
	// Create and insert a Site object
	UID uid = dao()->allocIdRangeForSites(1u);

	Site site;
	site.id = ++uid;
	site.name = "IRFU";
	dao()->insert(site);

	// Retrieve the site by name
	SitePtr sitePtr = dao()->searchSite("IRFU");
	ensure(sitePtr.get() != 0);
	ensure_equals(sitePtr->name, "IRFU");
}

/**
 * Test insertion of batch of sites.
 */
template<> template<>
void TestObject::test<14>()
{
	// Create a batch of 3 Site objects.
	SiteBatch batch;
	Site site;

	UID uid = dao()->allocIdRangeForSites(3u);

	site.id = ++uid;
	site.name = "CERN";
	batch.append(site);

	site.id = ++uid;
	site.name = "IN2P3";
	batch.append(site);

	site.id = ++uid;
	site.name = "GANIL";
	batch.append(site);

	// Insert the 3 previous sites.
	dao()->insert(batch);

	// Retrieve the inserted sites.
	SitePtrSetPtr sites = dao()->getSites();
}

/**
 * Test Run insertion.
 */
template<> template<>
void TestObject::test<15>()
{
	// Create and insert a Run object
	UID uid = dao()->allocIdRangeForRuns(1u);

	Run run;
	run.id = ++uid;
	SitePtr sitePtr = dao()->searchSite("IRFU");
	run.fkSite = sitePtr->getId();
	dao()->insert(run);

	// Retrieve the run
	RunPtrSetPtr runSetPtr = sitePtr->getRuns();
	ensure(runSetPtr.get() != 0);
}

/**
 * Test Run creation (method in Site class).
 */
template<> template<>
void TestObject::test<16>()
{

	SitePtr sitePtr = dao()->searchSite("GANIL");
	RunPtr run = sitePtr->createNewRun();
	dao()->insert(*run);

	// Retrieve the run
	RunPtrSetPtr runs = sitePtr->getRuns();
	ensure(runs.get() != 0);
}

/**
 * Test RunCondition insertion.
 */
template<> template<>
void TestObject::test<17>()
{

	// Retrieve run
	RunPtr run = dao()->searchRun("GANIL", 0);
	ensure(run.get() != 0);

	// Retrieve configuration version
	ConfigVersionPtr version = dao()->searchConfig("config1")->getLastVersion();
	ensure(version.get() != 0);

	// Create and insert run condition
	RunConditionPtr condition = run->createNewCondition(*version);
	dao()->insert(*condition);

	// Retrieve the run conditions
	RunConditionPtrSetPtr conditions = run->getConditions();
	ensure(conditions.get() != 0);

	// Check if we retrieved the right number of run conditions.
	ensure_equals(conditions->size(), 1u);

	// Retrieve configuration versions directly
	ConfigVersionPtrSetPtr versions = run->getVersions();
	ensure(versions.get() != 0);
	ensure_equals(versions->size(), 1u);
}

/**
 *
 */
//template<> template<>
//void TestObject::test<20>()
//{
//	dao()->dropAllTables();
//}

