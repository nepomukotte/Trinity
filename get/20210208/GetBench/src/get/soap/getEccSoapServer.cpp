// $Id: getEccSoapServer.cpp 1541 2016-03-14 12:49:18Z psizun $
#include "utl/Logging.h"
#include "mdaq/utl/CmdLineArgs.h"
#include "mdaq/utl/Server.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"
#include "get/soap/EccSoapServer.h"
#include "ecceccObject.h"
#include "ecc.nsmap"
#include ICEUTIL_HEADER(Thread.h)

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;

#include <algorithm>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <string>
#include <vector>

struct Namespace namespaces[] =
{   // {"ns-prefix", "ns-name"}
   {"SOAP-ENV", "http://schemas.xmlsoap.org/soap/envelope/", NULL, NULL}, // MUST be first
   {"SOAP-ENC", "http://schemas.xmlsoap.org/soap/encoding/", NULL, NULL}, // MUST be second
   {"xsi", "http://www.w3.org/2001/XMLSchema-instance", NULL, NULL}, // MUST be third
   {"xsd", "http://www.w3.org/2001/XMLSchema", NULL, NULL}, // Required for XML Schema types
   {"ecc", "urn:http://localhost:8083/", NULL, NULL}, // The namespace URI of the service operations
   {NULL, NULL, NULL, NULL} // end of table
};

static bool askedForShutdown = false;

get::soap::EccServer getEccServer;
//_________________________________________________________________________________________________
int ecc__ServerExit(struct soap *, ecc__Response & response)
{
	response.ErrorCode = 0;
	response.ErrorMessage = "";

	LOG_INFO() << "Shutting down ECC server...";

	askedForShutdown = true;

	return SOAP_OK;
}
//_________________________________________________________________________________________________
void sigint_handler(int)
{
	LOG_INFO() << "Shutting down ECC server...";
	askedForShutdown = true;
}
//_________________________________________________________________________________________________
class ServiceThread : public IceUtil::Thread
{
public:
	ServiceThread(struct soap* s) : s_(s) {}
	virtual ~ServiceThread()
	{
		soap_free(s_);    // Detach and deallocate context
	}
	// Overrides method from ::Ice::Thread
	virtual void run()
	{
		ecc_serve(s_);    // Process request
		soap_destroy(s_); // Clean up deserialized class instances
		soap_end(s_);     // Clean up deserialized data (except class instances) and temporary data
		soap_done(s_);    // Detach soap struct
		//soap_free(s_);    // Detach and deallocate context
	}
private:
	struct soap* s_;
};
typedef ::IceUtil::Handle < ServiceThread > ServiceThreadPtr;
//_________________________________________________________________________________________________
int main(int argc, char * argv[])
{
	try
	{
		// Default parameters
		int port=8083;
		uint16_t alarmPort;
		std::string repoUrl("."), repoUser, repoPwd;

		// Define options
		po::options_description allowed_options("Options");
		allowed_options.add_options()
			("help,h", "help message")
			("port,p", po::value<int>(&port), "ECC server port")
			("config-repo-url,h",   po::value<std::string>(&repoUrl)->default_value("."), "URL of the config. directory or database")
			("config-repo-user,u",   po::value<std::string>(&repoUser), "Username for database authentication")
			("config-repo-pwd,w",   po::value<std::string>(&repoPwd), "Password for database authentication")
			("alarm-port,a",   po::value<uint16_t>(&alarmPort)->default_value(get::Default::alarmLoggerPortNum), "Port for alarm notifications");
		po::positional_options_description p;
		p.add("address", 1);

		// Parse options
		po::variables_map vm;
		try
		{
			po::parsed_options parsed = po::command_line_parser(argc, argv).options(allowed_options).positional(p).allow_unregistered().run();
			po::store(parsed, vm);
			po::notify(vm);
		}
		catch (const po::error & e)
		{
			LOG_FATAL() << e.what();
			LOG_INFO() << "Usage: " << argv[0] << " [OPTIONS...]";
			LOG_INFO() << allowed_options << '\n';
			return EXIT_FAILURE;
		}

		// Process options
		if (vm.count("help"))
		{
			LOG_INFO() << "Usage: " << argv[0] << " [OPTIONS...]";
			LOG_INFO() << allowed_options << '\n';
			return EXIT_FAILURE;
		}
		getEccServer.setLoggingBackend();

		// Get port from environment variable ECC_PORT
		if (not vm.count("port"))
		{
			char* env_port = getenv("ECC_PORT");
			if (env_port != NULL)
			{
				port = atoi(env_port);
			}
		}

		// Disable IPv6 support
		// With Ice 3.5, support for IPv6 is enabled by default.
		// On sedipcc212 (Ubuntu 14.04), it causes "Address family is not supported by protocol" socket exceptions when using wildcard address.
		mdaq::utl::CmdLineArgs args = mdaq::utl::CmdLineArgs::build(argc, argv);
		args.push_back("--Ice.IPv6=0");

		//args.push_back("--Ice.Trace.ThreadPool=1");

		// Set ICE initialization parameters
		::mdaq::utl::Server::setInitData(args);

		// Setup configuration manager
		getEccServer.setConfigManager(repoUrl, repoUser, repoPwd);

		// Create alarm logger
		getEccServer.createAlarmLogger(alarmPort);

		eccService soap;
		soap.send_timeout = 60; // 60 seconds
		soap.recv_timeout = 60; // 60 seconds
		soap.accept_timeout = 3600; // server stops after 1 hour of inactivity
		soap.max_keep_alive = 100; // max keep-alive sequence

		// Bind
		LOG_INFO() << "Starting Electronics Control Core SOAP server on port " << port;
		soap.bind_flags |= SO_REUSEADDR;
		const int BACKLOG = 100; // Max. request backlog

		SOAP_SOCKET m = soap.bind(NULL, port, BACKLOG); // Returns master socket, backlog = max. queue size for requests
		if (!soap_valid_socket(m))
		{
			LOG_FATAL() << "Socket connection failed";
			return EXIT_FAILURE;
		}
		//LOG_DEBUG() << "Socket connection successful " << m;

		// Set handler for SIGINT signal to properly close server on Ctrl-C
		signal(SIGINT, sigint_handler);

		struct soap *tsoap;
		// Until server shutdown is requested
		while (not askedForShutdown)
		{
			// Wait for new client
			soap.accept_timeout = 5;
			SOAP_SOCKET s = soap_accept(&soap); // Returns slave socket
			if (!soap_valid_socket(s))
			{
				if (soap.errnum)
				{
					soap_print_fault(&soap, stderr);
				}
				continue; // Retry
			}
			//LOG_DEBUG() << "Accepted socket " << s << " connection from IP "
			//	<< ((service.ip >> 24) & 0xFF) << '.' <<  ((service.ip >> 16) & 0xFF) << '.'
			//	<< ((service.ip >> 8)  & 0xFF) << '.' << (service.ip & 0xFF);


			tsoap = soap_copy(&soap); // Make a safe copy
			tsoap->socket = s; // new socket fd

			// Handle request in new thread
			ServiceThreadPtr t = new ServiceThread(tsoap);
			t->start().detach();
		}
		// Detach and deallocate context
		soap_done(&soap);

		LOG_DEBUG() << "Bye";

		// Switch off log4cxx logging backend (see issue LOGCXX-322)
		::utl::BackendLogger::setBackend(::utl::LoggingBackendPtr(new ::mdaq::utl::ConsoleLoggingBackend));

		// Cleanup ICE
		::mdaq::utl::Server::ic()->destroy();

		return EXIT_SUCCESS;
	}
	catch (const std::exception & e)
	{
		LOG_FATAL() << e.what();
	}

	return EXIT_FAILURE;
}
//_________________________________________________________________________________________________
int ecc__GetConfigIDs(struct soap *, ecc__ResponseText & respText)
{
	return getEccServer.getConfigIDs(respText);
}
//_________________________________________________________________________________________________
int ecc__Describe(struct soap *, std::string xmlConfigId, std::string options, ecc__ResponseText & response)
{
	return getEccServer.describe(xmlConfigId, options, response);
}
//_________________________________________________________________________________________________
int ecc__Prepare(struct soap *, std::string xmlConfigId, std::string options, ecc__ResponseText & response)
{
	return getEccServer.prepare(xmlConfigId, options, response);
}
//_________________________________________________________________________________________________
int ecc__Undo(struct soap *, std::string , std::string, ecc__ResponseText & response)
{
	return getEccServer.undo(response);
}
//_________________________________________________________________________________________________
int ecc__Configure(struct soap *, std::string xmlConfigId, std::string dataLinks, ecc__ResponseText & response)
{
	return getEccServer.configure(xmlConfigId, dataLinks, response);
}
//_________________________________________________________________________________________________
int ecc__Start(struct soap *, std::string , std::string, ecc__ResponseText & response)
{
	return getEccServer.start(response);
}
//_________________________________________________________________________________________________
int ecc__Stop(struct soap *, std::string , std::string, ecc__ResponseText & response)
{
	return getEccServer.stop(response);
}
//_________________________________________________________________________________________________
int ecc__Pause(struct soap *, std::string , std::string, ecc__ResponseText & response)
{
	return getEccServer.pause(response);
}
//_________________________________________________________________________________________________
int ecc__Resume(struct soap *, std::string , std::string, ecc__ResponseText & response)
{
	return getEccServer.resume(response);
}
//_________________________________________________________________________________________________
int ecc__Breakup(struct soap *, std::string , std::string, ecc__ResponseText & response)
{
	return getEccServer.breakup(response);
}
//_________________________________________________________________________________________________
int ecc__GetState(struct soap *, ecc__ResponseState & response)
{
	return getEccServer.getState(response);
}
//_________________________________________________________________________________________________
int ecc__Special(struct soap *, std::string command, std::string options, ecc__ResponseText & respText)
{
	return getEccServer.specialCmd(command, options, respText);
}
//_________________________________________________________________________________________________
