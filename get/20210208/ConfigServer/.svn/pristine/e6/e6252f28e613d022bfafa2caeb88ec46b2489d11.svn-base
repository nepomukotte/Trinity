/*
 * ConfigServer.cpp
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

BOOST_AUTO_TEST_CASE( test_loadDoc )
{
	cout << "Testing ConfigServer:" << endl;
	RDocumentPrx mydoc = serv->loadDoc("test/ConfigTest.xml");
	RCConfigPrx conf = mydoc->getConfig();
	BOOST_CHECK_EQUAL(conf->getName(""), "Test");
}

BOOST_AUTO_TEST_CASE( test_newDoc )
{
	RDocumentPrx mydoc = serv->newDoc();
	RCConfigPrx conf = mydoc->getConfig();
	BOOST_CHECK_EQUAL(conf->getName(""), "unnamed");
}

BOOST_AUTO_TEST_SUITE_END()


//____________________________________________________________________________//
