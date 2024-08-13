/* =====================================================================================================================
 * getEccServer.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * getEccServer functions
 * Created on: Jan 11, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of the Mordicus Real-Time Software
 * Development Framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#include "get/GetEccImpl.h"
#include "get/mt/AlarmLogger.h"
#include "mdaq/DefaultPortNums.h"
#include "mdaq/utl/Server.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"
#ifdef HAVE_ICESTORM
#include "mdaq/utl/IceStormLoggingBackend.h"
#endif
#ifdef HAVE_LOG4CXX
#include "mdaq/utl/Log4cxxLoggingBackend.h"
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/defaultconfigurator.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/logmanager.h>
#include <log4cxx/patternlayout.h>
#endif
#include "CCfg/Io/Url.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <sstream>

using ::mdaq::utl::Server;
using ::mdaq::utl::CmdLineArgs;
using ::get::GetEccImpl;

static Server* getEccServer;

enum LoggingBackendType { ConsoleLogger, IceStormLogger, Log4cxxLogger };

std::istream& operator>>(std::istream & in, LoggingBackendType & type)
{
    std::string token;
    in >> token;
    if (token == "LOG4CXX")
    	type = Log4cxxLogger;
    else if (token == "ICESTORM")
    	type = IceStormLogger;
    else if (token == "CONSOLE")
    	type = ConsoleLogger;
    else
    	// throw po::validation_error(po::validation_error::invalid_option_value); // no such constructor in Boost 1.41
    	throw po::error("invalid logging backend type");
    return in;
}

int main(int argc, char* argv[])
{
	try
	{
		// Define options
		std::string repoUrl("."), repoUser, repoPwd, serverAddress, alarmsAddress;
		LoggingBackendType loggerType = ConsoleLogger;
		po::options_description allowed_options("Options");
		allowed_options.add_options()("help", "help message")
			 ("address,a", po::value<std::string>(&serverAddress), "ECC server host address and port")
    	     ("config-repo-url,h",   po::value<std::string>(&repoUrl)->default_value("."), "URL of the config. directory or database")
    	     ("config-repo-user,u",   po::value<std::string>(&repoUser), "Username for database authentication")
    	     ("config-repo-pwd,p",   po::value<std::string>(&repoPwd), "Password for database authentication")
    	     ("logger,l", po::value< LoggingBackendType >(&loggerType), "Logging backend (CONSOLE, LOG4CXX, ICESTORM)")
    	     ("alarm-address,r",   po::value<std::string>(&alarmsAddress), "Address and port for alarm notifications");
		po::positional_options_description p;
		p.add("address", 1);

		po::variables_map vm;
		CmdLineArgs args;
		try
		{
			// Parse options
			po::parsed_options parsed = po::command_line_parser(argc, argv).options(allowed_options).positional(p).allow_unregistered().run();
			po::store(parsed, vm);
			po::notify(vm);
			args = po::collect_unrecognized(parsed.options, po::include_positional);
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

		// Disable IPv6 support
		// With Ice 3.5, support for IPv6 is enabled by default.
		// On sedipcc212 (Ubuntu 14.04), it causes "Address family is not supported by protocol" socket exceptions when using wildcard address.
		args.push_back("--Ice.IPv6=0");

		// Set ICE communicator initialization data
		::mdaq::utl::Server::setInitData(args);

		// Setup logging backend
		utl::BackendLogger::setBackend(utl::LoggingBackendPtr(new mdaq::utl::ConsoleLoggingBackend()));
		if (IceStormLogger == loggerType)
		{
#ifdef HAVE_ICESTORM
			utl::BackendLogger::setBackend(utl::LoggingBackendPtr(new mdaq::utl::IceStormLoggingBackend("GetEccLogger")));
#endif
		}
		else if (Log4cxxLogger == loggerType)
		{
#ifdef HAVE_LOG4CXX
		// Initiate default initialization procedure.
		log4cxx::DefaultConfigurator::configure(log4cxx::LogManager::getLoggerRepository());
		// Check if GET logger has been defined by log4j property file
		if (log4cxx::LogManager::getLogger("GET")->getAllAppenders().empty())
		{
			// Configure custom console logger
			log4cxx::LoggerPtr logger = log4cxx::LogManager::getLogger("GET");
			log4cxx::ConsoleAppenderPtr appender(new log4cxx::ConsoleAppender());
			log4cxx::LayoutPtr layout(new log4cxx::PatternLayout("%d{(%H:%M:%S - %d/%m/%Y)} %p - %m%n"));
			appender->setLayout(layout);
			log4cxx::helpers::Pool pool;
			appender->activateOptions(pool);
			log4cxx::Logger::getRootLogger()->addAppender(appender);
			log4cxx::LogManager::getLoggerRepository()->setConfigured(true);
		}

		// Setup logging backend
		utl::BackendLogger::setBackend(utl::LoggingBackendPtr(new mdaq::utl::Log4cxxLoggingBackend("GET.ECC")));
#endif
		}

		LOG_INFO() << "Starting Electronics Control Core server...";

		// Create server
		get::GetEccImplPtr getEccImpl;
		mdaq::utl::net::SocketAddress sockAddr(INADDR_ANY, ::mdaq::Default::eccServerCtrlPortNum);
		if (vm.count("address"))
		{
			sockAddr.fromString(serverAddress);
			if (sockAddr.isValid())
			{
				if (sockAddr.port() == 0)
				{
					sockAddr.port() = ::mdaq::Default::eccServerCtrlPortNum;
				}
				getEccImpl = new get::GetEccImpl();
				getEccServer = &Server::create(sockAddr).addServant("Ecc", getEccImpl);
			}
			else
			{
				LOG_FATAL() << "Invalid IP address / port: " << sockAddr;
				return EXIT_FAILURE;
			}
		}
		else
		{
			getEccImpl = new get::GetEccImpl();
			getEccServer = &Server::create(sockAddr).addServant("Ecc", getEccImpl);
		}

		// Setup configuration manager
		if (vm.count("config-repo-url"))
		{
			CCfg::Io::Url repoFullUrl(repoUrl);
			if (not repoFullUrl.isFile())
			{
				if (not repoUser.empty())
				{
					repoFullUrl.setUser(repoUser);
				}
				if (not repoPwd.empty())
				{
					repoFullUrl.setPassword(repoPwd);
				}
			}
			getEccImpl->setConfigManager(repoFullUrl.format());
		}

		// Create alarm logger
		::utl::net::SocketAddress alarmLoggerAddr(sockAddr);
		if (vm.count("alarm-address"))
		{
			alarmLoggerAddr.fromString(alarmsAddress);
		}
		else if (alarmLoggerAddr.ipAddress().in_addr() != INADDR_ANY)
		{
			// Ensure alarm logger is available from all IP addresses
			alarmLoggerAddr.ipAddress().in_addr() = INADDR_ANY;
			alarmLoggerAddr.port() = get::Default::alarmLoggerPortNum;
		}
		getEccImpl->createAlarmLogger(alarmLoggerAddr);

		// Start server thread
		getEccServer->start();
		// Wait for shutdown
		getEccServer->ic()->waitForShutdown();

		return EXIT_SUCCESS;
	}
	catch (const std::exception& e) { LOG_ERROR() << e.what(); }
	catch (const std::string&  msg) { LOG_ERROR() << msg; }
	catch (const char*         msg) { LOG_ERROR() << msg; }
	catch (...)                     { LOG_ERROR() << "Unknown exception"; }

	if (getEccServer)
	{
		getEccServer->adapter()->destroy();
	}

	return EXIT_FAILURE;
}
