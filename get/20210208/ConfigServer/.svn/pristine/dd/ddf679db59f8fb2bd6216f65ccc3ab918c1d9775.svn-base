/*
 * RDocumentTest.cpp
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


BOOST_AUTO_TEST_CASE( test_getConfig )
{
	cout << "Testing RDocument:" << endl;
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	BOOST_CHECK_EQUAL(conf->getName(""), "Test");
}

BOOST_AUTO_TEST_CASE( test_saveToFile )
{
	RDocumentPrx mydoc = serv->newDoc();
	RCConfigPrx conf = mydoc->getConfig();
	conf->setInt("a", "42");
	mydoc->saveToFile("test/SavedConfig.xml");
	mydoc = serv->loadDoc("test/SavedConfig.xml");
	conf = mydoc->getConfig();
	BOOST_CHECK_EQUAL(conf->getInt("a"), 42);
}

BOOST_AUTO_TEST_CASE( test_saveToUrl )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->removeNode("Numbers");
	mydoc->saveToUrl("test/SavedConfig1.xml"); //Could be used with a Url to a Database
	mydoc = serv->loadDoc("test/SavedConfig1.xml");
	BOOST_CHECK_EQUAL(conf->testNode("Numbers"), 0);
}

BOOST_AUTO_TEST_CASE( test_save )
{
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	conf->removeNode("Numbers");
	mydoc->saveToFile("test/SavedConfig2.xml");
	mydoc = serv->loadDoc("test/SavedConfig2.xml");
	conf = mydoc->getConfig();
	BOOST_CHECK_EQUAL(conf->testNode("Numbers"), 0);
}

BOOST_AUTO_TEST_SUITE_END()


//____________________________________________________________________________//
