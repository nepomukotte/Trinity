#include "mdaq/EccCmd.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"
#include "mdaq/utl/Server.h"
#include "mdaq/utl/CmdLineArgs.h"
#ifdef HAVE_ICESTORM
#include "mdaq/utl/IceStormLoggingSubscriber.h"
#endif
#include "mdaq/DefaultPortNums.h"

using namespace mdaq::utl;

int main(int argc, char** argv)
{
	// Set ICE communicator initialization data
	mdaq::utl::Server::setInitData(argc, argv);

	// Setup logging backend
	utl::BackendLogger::setBackend(utl::LoggingBackendPtr(new mdaq::utl::ConsoleLoggingBackend()));

	// Subscribe to ECC server remote logger
#ifdef HAVE_ICESTORM
	IceStormLoggingSubscriber remoteLogger("MDaqEccLogger");
#endif

	// Run command line client
	mdaq::EccCmd<>::Engine eccClientCmdLine;
	int status = eccClientCmdLine.main(mdaq::utl::Server::initData());

	// Destroy communicator if any
	try
	{
		mdaq::utl::Server::ic()->destroy();
	}
	catch (const Ice::Exception & e) {}
	return status;
}

