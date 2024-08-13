#include "get/GetTest.h"
#include "get/WorkspaceWizard.h"
#include "get/WorkspaceManager.h"
#include "get/daq/DataRouterWizard.h"
#include "get/DefaultPortNums.h"
#include "get/ControllerOptions.h"
#include "mdaq/utl/CmdLineArgs.h"
#include "mdaq/utl/Server.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"
#ifdef HAVE_ICESTORM
#include "mdaq/utl/IceStormLoggingSubscriber.h"
#endif
#include <mfm/FrameDictionary.h>
using namespace get;
using mdaq::utl::CmdLineArgs;

#include <QApplication>
#include <QSettings>
#include <QIcon>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <cstdlib>
#include <cstring>
//__________________________________________________________________________________________________
bool load_format(const std::string & userFormats, const std::string & defaultFormats="CoboFormats.xcfg")
{
	LOG_INFO() << "Loading frame formats...";
	if (not userFormats.empty())
	{
		try
		{
			LOG_DEBUG() << "Loading formats from file '" << userFormats << '\'';
			mfm::FrameDictionary::instance().addFormats(userFormats);
		}
		catch (const CCfg::Exception & e)
		{
			LOG_FATAL() << e.what();
			return EXIT_FAILURE;
		}
	}
	else
	{
		try
		{
			LOG_DEBUG() << "Loading formats from file '" << defaultFormats << '\'';
			mfm::FrameDictionary::instance().addFormats(defaultFormats);
		}
		catch (const CCfg::Exception & e)
		{
			LOG_WARN() << e.what();
			try
			{
				const std::string officialFormats = QString("%1/../share/get-bench/format/CoboFormats.xcfg").arg(QCoreApplication::applicationDirPath()).toStdString();
				LOG_DEBUG() << "Loading formats from file '" << officialFormats << '\'';
				mfm::FrameDictionary::instance().addFormats(officialFormats);
			}
			catch (const CCfg::Exception & e)
			{
				LOG_FATAL() << e.what() << "\nCould not load frame formats!";
				return EXIT_FAILURE;
			}
		}
	}
	return EXIT_SUCCESS;
}
//__________________________________________________________________________________________________
int main(int argc, char *argv[])
{
	std::string formatsFile;
	uint16_t alarmPort;

	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("format,f", po::value<std::string>(&formatsFile), "overrides format description file")
		("alarm-port,a",   po::value<uint16_t>(&alarmPort), "Port for alarm notifications")
	;
	po::positional_options_description p;

	// Parse options
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		LOG_INFO() << "Usage: " << argv[0] << " [OPTIONS...]";
		LOG_INFO() << desc;
		return EXIT_FAILURE;
	}

	// Disable IPv6 support
	// With Ice 3.5, support for IPv6 is enabled by default.
	// On sedipcc212 (Ubuntu 14.04), it causes "Address family is not supported by protocol" socket exceptions when using wildcard address.
	CmdLineArgs args = CmdLineArgs::build(argc, argv);
	args.push_back("--Ice.IPv6=0");

    // Set ICE initialization parameters
    ::mdaq::utl::Server::setInitData(args);

    // Set logging backend
	::utl::BackendLogger::setBackend(::utl::LoggingBackendPtr(new ::mdaq::utl::ConsoleLoggingBackend));

	// Subscribe to ECC server remote logger
#ifdef HAVE_ICESTORM
	::mdaq::utl::IceStormLoggingSubscriber remoteLogger("GetEccLogger");
#endif

	// Create Qt application
    QApplication app(argc, argv);
    // Load resources
    Q_INIT_RESOURCE(ConfigWidgets);
    Q_INIT_RESOURCE(GetController);
    QApplication::setApplicationName(QObject::tr("GetController"));
    QApplication::setWindowIcon(QPixmap(":/get/logo-highres.png"));
    QApplication::setOrganizationName("actar-get");
    QApplication::setOrganizationDomain("www-actar-get.cea.fr");

    // Load data formats
    if (EXIT_SUCCESS != load_format(formatsFile))
    	return EXIT_FAILURE;

	// Theme
#if QT_VERSION >= 0x040600
	if (not QIcon::hasThemeIcon("application-exit"))
	{
		// If there is no default working icon theme (Windows, Mac OS X)
		// then we provide fallback theme
		QIcon::setThemeName("oxygen");
	}
#endif

    // Select workspace
	QString workspacePath = WorkspaceWizard::selectWorkspace(0, false);
	if (workspacePath.isEmpty())
		return EXIT_SUCCESS;
	if (not WorkspaceManager::instance().loadWorkspace(workspacePath.toStdString()))
		return EXIT_FAILURE;

	// Create main GUI window
    GetTest mainWindow;
    mainWindow.show();
    mainWindow.raise();

    // Setup data router
    get::daq::DataRouterWizard wizard;
    if (wizard.exec() != QDialog::Accepted)
    	return EXIT_SUCCESS;
    mainWindow.createDataRouter();

    // Create alarm logger
    utl::net::SocketAddress alarmAddr = ControllerOptions::instance().alarmLoggerEndpoint();
    if (vm.count("alarm-port"))
    {
    	alarmAddr.port() = alarmPort;
    }
    mainWindow.createAlarmLogger(alarmAddr);

    // Enter event loop
    int status = app.exec();

	// Destroy communicator if any
	try
	{
		mdaq::utl::Server::ic()->destroy();
	}
	catch (const Ice::Exception &) {}

	return status;
}
