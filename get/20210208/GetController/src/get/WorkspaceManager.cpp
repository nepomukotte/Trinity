/**
 * $Id: WorkspaceManager.cpp 1200 2014-03-20 17:08:56Z psizun $
 * @file WorkspaceManager.cpp
 * @date 19 avr. 2012
 * @author sizun
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * -----------------------------------------------------------------------------
 */

#include "get/WorkspaceManager.h"
#include "get/GetControllerException.h"
#include "mdaq/DefaultPortNums.h"
#include "mdaq/utl/net/SocketAddress.h"
#include "utl/Logging.h"
#include <CCfg/Scalar.h>
#include <CCfg/Range.h>
#include <CCfg/Enum.h>
#include <CCfg/View/ClassIterator.h>
#include <CCfg/Io/Url.h>
using namespace CCfg;

#include <QString>
#include <QDir>
#include <QFileInfo>
#include <vector>
#include <sstream>
#include <stdexcept>

namespace get {
//_____________________________________________________________________________
/**
 * Allowed test types:
 *   - Default (formerly BaseLine): GetController is passive during run.
 *   - StdAcquisition: depreciated since embedded software provides periodic pulser.
 *   - Acq2pMode: idem in two proton decay mode.
 *   - PeakingTime: generates N events for each of the 16 peaking time values.
 */
const char* WorkspaceManager::ALLOWED_TEST_TYPES[] = {"Acq2pMode", "Default", "PeakingTime", "StdAcquisition"};
//_____________________________________________________________________________
WorkspaceManager::WorkspaceManager()
{
}
//_____________________________________________________________________________
WorkspaceManager::~WorkspaceManager()
{
}
//_____________________________________________________________________________
WorkspaceManager & WorkspaceManager::instance()
{
	static WorkspaceManager testManager;
	return testManager;
}
//_____________________________________________________________________________
CCfg::Io::Document & WorkspaceManager::document() const
{
	if (not testsDoc_.get())
	{
		throw GetControllerException("No workspace has been loaded");
	}
	return *testsDoc_.get();
}
//_____________________________________________________________________________
CCfg::CConfig WorkspaceManager::workspace() const
{
	if (not testsDoc_.get())
	{
		throw GetControllerException("No workspace has been loaded");
	}
	return CConfig(testsDoc_->getConfig());
}
//_____________________________________________________________________________
CCfg::CConfig WorkspaceManager::workspace(const std::string & name, const std::string & index) const
{
	return workspace()(name, index);
}
//_____________________________________________________________________________
CCfg::CConfig WorkspaceManager::defaultTest()
{
	return workspace("Test");
}
//_____________________________________________________________________________
void WorkspaceManager::createWorkspace(const std::string & path)
{
	LOG_DEBUG() << "Creating GetController workspace \"" << path << '"';

	// Create config from scratch
	std::auto_ptr< CCfg::Object > newWorkspace(new CCfg::Object("Workspace"));

	// Create default test
	CCfg::CConfig workspaceCfg(newWorkspace.get());
	createDefaultTest(workspaceCfg);

	// Save document
	CCfg::Io::Document doc;
	doc.setConfig(newWorkspace);
	doc.saveToFile(path, false);
}
//_____________________________________________________________________________
void WorkspaceManager::createDefaultTest(CCfg::CConfig & config)
{
	try {
		// Create default test
		Attribute & workspace = config.find();
		workspace.requireObject("Test", "*");
		CConfig defaultTest = config("Test", "*");
		defaultTest("Type").setStringValue("StdAcquisition");
		CCfg::Enum allowedTypes;
		size_t const typeCount = sizeof(ALLOWED_TEST_TYPES)/sizeof(*ALLOWED_TEST_TYPES);
		for (size_t i=0; i < typeCount; ++i)
		{
			allowedTypes.insert(ALLOWED_TEST_TYPES[i]);
		}
		defaultTest("Type").setRange(allowedTypes);
		defaultTest("Setup", "*").create(OBJECT_INST);
		defaultTest("Setup", "Hardware").create(OBJECT_INST);
		defaultTest("Setup", "Hardware")("Path").setStringValue("");
		defaultTest("Setup", "Conditions").create(OBJECT_INST);
		defaultTest("Setup", "Conditions")("Path").setStringValue("");
		defaultTest("Setup", "Data")("Path").setStringValue("/tmp");

		// Setup default network
		defaultTest.find().requireObject("Setup", "Network").requireObject("IpNode", "*");
		CConfig defaultSocketAddress = defaultTest("Setup", "Network")("IpNode", "*");
		mdaq::utl::net::SocketAddress("192.168.0.0", mdaq::Default::hwServerCtrlPortNum).toCConfig(defaultSocketAddress);
		defaultSocketAddress("addressField", "*").setIntValue(255);
		defaultSocketAddress("addressField", "*").setRange(Range::createFromString("[0, 255]"));
		defaultSocketAddress("port").setRange(Range::createFromString("[0, 65535]"));

		// Set default ports
		CConfig network = defaultTest("Setup", "Network");
		network("IpNode", "ECC")("port").setIntValue(mdaq::Default::eccServerCtrlPortNum);
		network("IpNode", "DAQ")("port").setIntValue(mdaq::Default::dataRouterCtrlPortNum);
		network("IpNode", "TARGET").create(OBJECT_INST);
	}
	catch (const CCfg::Exception & e)
	{
		LOG_ERROR() << "Could not create default test: " << e.what();
	}
}
//_____________________________________________________________________________
void WorkspaceManager::saveTest(const TestSummary & summary)
{
	summary.toCConfig(workspace());
	save();
}
//_____________________________________________________________________________
bool WorkspaceManager::loadWorkspace(const std::string & path)
{
	closeCurrentTest();

	LOG_INFO() << "Loading workspace \"" << path << '"';
	testsDoc_.reset(new CCfg::Io::Document());
	try
	{
		testsDoc_->loadFromFile(path);
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
		return false;
	}

	return loadTests();
}
//_____________________________________________________________________________
/**
 * Parses document and stores the name of all tests.
 */
bool WorkspaceManager::loadTests()
{
	LOG_DEBUG() << "Searching for tests in workspace ...";
	// Reset list of tests
	testNames_.clear();

	// Check document is open
	if (not testsDoc_.get())
		return false;

	// Find all concrete tests with class path Workspace[*].Test[*]
	try
	{
		View::ClassIterator testIter = workspace("Test").getClassInstances();
		while (testIter.hasNext())
		{
			// Do not forget to call next before accessing the first element.
			testIter.next();
			// Test whether object is concrete
			if (not (testIter->toObject().isDefaultConfig()))
			{
				testNames_.insert(testIter->getIndex());
			}
		}
		LOG_DEBUG() << "Found " << testNames_.size() << " test(s)";
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		LOG_FATAL() << "It really does not look like a valid workspace file. " << e.what();
		testsDoc_.reset();
		return false;
	}

	// Set current test
	if (testNames_.find(currentTestName_) == testNames_.end())
	{
		if (testNames_.size() > 0)
			currentTestName_ = *testNames_.begin();
		else
			currentTestName_.clear();
	}

	return true;
}
//_____________________________________________________________________________
size_t WorkspaceManager::testCount() const
{
	return testNames().size();
}
//_____________________________________________________________________________
/**
 * Parses document and stores the name of all tests.
 */
void WorkspaceManager::printTests(std::ostream & stream)
{
	std::set< std::string >::iterator testNameIter;
	for (testNameIter = testNames_.begin(); testNameIter != testNames_.end(); testNameIter++)
	{
		TestSummary descr(workspace("Test", *testNameIter));
		descr.print(stream);
	}
}
//_____________________________________________________________________________
std::auto_ptr< ::get::TestSummary > WorkspaceManager::currentTest() const
{
	return readTest(currentTestName());
}
//_____________________________________________________________________________
void WorkspaceManager::setCurrentTest(const std::string & name)
{
	if (testNames_.find(name) != testNames_.end())
	{
		currentTestName_ = name;
	}
}
//_____________________________________________________________________________
void WorkspaceManager::closeCurrentTest()
{
	currentTestName_.clear();
}
//_____________________________________________________________________________
/**
 * Reads CompoundConfig document ands returns test summary.
 *
 */
//_____________________________________________________________________________
std::auto_ptr< ::get::TestSummary > WorkspaceManager::readTest(const std::string & name) const
{
	std::auto_ptr< get::TestSummary > summaryPtr(new TestSummary(workspace("Test", name)));
	return summaryPtr;
}
//_____________________________________________________________________________
/**
 * Saves current document.
 */
void WorkspaceManager::save()
{
	if (not testsDoc_.get())
		return;
	testsDoc_->saveToFile(testsDoc_->getUrl(), false);
	// Reload tests
	loadTests();
}
//_____________________________________________________________________________
} /* namespace get */
