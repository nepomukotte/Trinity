/*
 * RCConfigTest.cpp
 *
 *  Created on: 7 juil. 2009
 *      Author: rvloic
 */

#include <boost/test/auto_unit_test.hpp>
#include <iostream>
#include ICE_HEADER(Initialize.h)
#include ICE_HEADER(Communicator.h)
#include "ConfigServer.h"

using namespace std;
using namespace CCfg;


//____________________________________________________________________________//

struct TestServer {
public:
    Ice::CommunicatorPtr ic;
    Ice::ObjectPrx base;
    ConfigServerPrx serv;
public:
	TestServer() {
		int dummy = 0;
		ic = Ice::initialize(dummy, 0);
		base = this->ic->stringToProxy("ConfigServer:default -p 10000");
		serv = ConfigServerPrx::checkedCast(base);
	}
	~TestServer()  { }
};

//____________________________________________________________________________//

BOOST_FIXTURE_TEST_SUITE( RDocument_suite, TestServer)
/////////////////////////

//Config Functions
BOOST_AUTO_TEST_CASE( test_append )
{
	cout << "Testing RCConfig:" << endl;
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->append("Numbers");
	conf->removeNode("");
	BOOST_CHECK_EQUAL(conf->testNode("Numbers"), 0);
}

BOOST_AUTO_TEST_CASE( test_subConf )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	RCConfigPrx conf2 = conf->subConf("Numbers");
	BOOST_CHECK_EQUAL(conf2->getName(""), "Numbers");
}

BOOST_AUTO_TEST_CASE( test_testNode )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	BOOST_CHECK_EQUAL(conf->testNode("Numbers"), 1);
}

BOOST_AUTO_TEST_CASE( test_createNode )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->createNode("NewNode");
	BOOST_CHECK_EQUAL(conf->testNode("NewNode"), 1);
}

BOOST_AUTO_TEST_CASE( test_removeNode )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->removeNode("Numbers");
	BOOST_CHECK_EQUAL(conf->testNode("Numbers"), 0);
}

BOOST_AUTO_TEST_CASE( test_enumNode )
{
	RAttributeList	nodeList;

	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->enumNode("", nodeList);
	BOOST_CHECK_EQUAL(nodeList[0], "Numbers");
}

//Get Functions
BOOST_AUTO_TEST_CASE( test_getCurrentPath )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->append("Numbers");
	BOOST_CHECK_EQUAL(conf->getCurrentPath(), "Numbers");
}

BOOST_AUTO_TEST_CASE( test_getInt )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->append("Numbers");
	BOOST_CHECK_EQUAL(conf->getInt("i"), 3);
}

BOOST_AUTO_TEST_CASE( test_getStr )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setStr("str", "yop");
	BOOST_CHECK_EQUAL(conf->getStr("str"), "yop");
}

BOOST_AUTO_TEST_CASE( test_getBool )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setBool("bool", "true");
	BOOST_CHECK_EQUAL(conf->getBool("bool"), true);
}

BOOST_AUTO_TEST_CASE( test_getReal )
{
	float temp = 4.2;
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setReal("real", "4.2");
	BOOST_CHECK_EQUAL(conf->getReal("real"), temp);
}

BOOST_AUTO_TEST_CASE( test_getName )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	BOOST_CHECK_EQUAL(conf->getName(""), "Test");
}

BOOST_AUTO_TEST_CASE( test_getUnit )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->append("Numbers");
	BOOST_CHECK_EQUAL(conf->getUnit("memSize"), "bytes");
}

BOOST_AUTO_TEST_CASE( test_getIndex )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	BOOST_CHECK_EQUAL(conf->getIndex("ind[0]"), "0");
}

BOOST_AUTO_TEST_CASE( test_getRange )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setRange("test", "[0, 100]");
	BOOST_CHECK_EQUAL(conf->getRange("test"), "[0, 100]");
}

BOOST_AUTO_TEST_CASE( test_getPreComment )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setPreComment("test", "PreCom");
	BOOST_CHECK_EQUAL(conf->getPreComment("test"), "PreCom");
}

BOOST_AUTO_TEST_CASE( test_getPostComment )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setPostComment("test", "PostCom");
	BOOST_CHECK_EQUAL(conf->getPostComment("test"), "PostCom");
}

//Set Functions
BOOST_AUTO_TEST_CASE( test_setCurrentPath )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setCurrentPath("Numbers");
	BOOST_CHECK_EQUAL(conf->getCurrentPath(), "Numbers");
}

BOOST_AUTO_TEST_CASE( test_setInt )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->append("Numbers");
	conf->setInt("j", "5");
	BOOST_CHECK_EQUAL(conf->getInt("j"), 5);
}

BOOST_AUTO_TEST_CASE( test_setStr )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setStr("str", "yop");
	BOOST_CHECK_EQUAL(conf->getStr("str"), "yop");
}

BOOST_AUTO_TEST_CASE( test_setBool )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setBool("bool", "true");
	BOOST_CHECK_EQUAL(conf->getBool("bool"), true);
}

BOOST_AUTO_TEST_CASE( test_setReal )
{
	float temp = 4.2;
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setReal("real", "4.2");
	BOOST_CHECK_EQUAL(conf->getReal("real"), temp);
}

BOOST_AUTO_TEST_CASE( test_setName )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setName("", "Name");
	BOOST_CHECK_EQUAL(conf->getName(""), "Name");
}

BOOST_AUTO_TEST_CASE( test_setUnit )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->append("Numbers");
	conf->setUnit("memSize", "meters");
	BOOST_CHECK_EQUAL(conf->getUnit("memSize"), "meters");
}

BOOST_AUTO_TEST_CASE( test_setIndex )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setIndex("ind[0]", "1");
	BOOST_CHECK_EQUAL(conf->getIndex("ind[1]"), "1");
}

BOOST_AUTO_TEST_CASE( test_setRange )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setRange("test", "[0, 100]");
	BOOST_CHECK_EQUAL(conf->getRange("test"), "[0, 100]");
}

BOOST_AUTO_TEST_CASE( test_setPreComment )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setPreComment("test", "PreCom");
	BOOST_CHECK_EQUAL(conf->getPreComment("test"), "PreCom");
}

BOOST_AUTO_TEST_CASE( test_setPostComment )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->setPostComment("test", "PostCom");
	BOOST_CHECK_EQUAL(conf->getPostComment("test"), "PostCom");
}

///////////////////////////
BOOST_AUTO_TEST_SUITE_END()


//____________________________________________________________________________//
