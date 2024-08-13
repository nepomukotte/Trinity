// $Id: cobo-frame-viewer.cpp 1667 2017-07-06 12:58:35Z psizun $
#include "get/CoBoFrameViewer.h"
#include "get/CoBoEvent.h"
#include "get/ChannelMask.h"
#include <CCfg/Exception.h>
#include <mfm/FrameDictionary.h>
#include <utl/Logging.h>
using namespace mfm;
using namespace get;

#include <QObject>
#include <QApplication>
#include <QString>
#include <QIcon>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
	int status = EXIT_SUCCESS;

	try
	{
		// List of options
		std::string filename;
		size_t frameIndex = 0;
		size_t agetIdx = 0u;
		std::string channelsString = "0-67";
		std::string fmtfile;
		size_t byteOffset = 0;

		// Define options
		po::options_description allowed_options("Options");
		allowed_options.add_options()("help,h", "help message")
    	     ("input-file,i",   po::value<std::string>(&filename), "filename")
			 ("offset,o",   po::value<size_t>(&byteOffset), "byte offset")
    	     ("frame,n",   po::value<size_t>(&frameIndex), "index of the frame within the file")
			 ("aget,a", po::value<size_t>(&agetIdx)->default_value(agetIdx), "AGET index")
			 ("channels,c", po::value<std::string>(&channelsString)->default_value(channelsString), "channels (ex: 1-5,7-63)")
			 ("format,f", po::value<std::string>(&fmtfile), "overrides format description file")
			 ("play,p", "start automatic scan of events")
			 ("3d", "enable 3D display mode");
		po::positional_options_description p;
		p.add("input-file", -1);

		// Parse options
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).
		          options(allowed_options).positional(p).run(), vm);
		po::notify(vm);
		// Process options
		if (vm.count("help"))
		{
			LOG_INFO() << "Usage: " << argv[0] << " -i frame.graw [OPTIONS...]";
			LOG_INFO() << allowed_options << '\n';
			status = EXIT_FAILURE;
			return status;
		}
		get::ChannelMask channels;
		if (vm.count("channels"))
		{
			channels = ChannelMask::fromString(channelsString);
		}

		std::ostringstream oss;
		std::ostream_iterator<int> osi(oss," ");
		std::copy(channels.begin(), channels.end(), osi);
		LOG_INFO() << "Selected channel(s): " << oss.str();


		QApplication app(argc, argv);
		Q_INIT_RESOURCE(CoBoFrameViewer);
		QApplication::setWindowIcon(QPixmap(":/get/logo.png"));
	    QApplication::setOrganizationName("actar-get");
	    QApplication::setOrganizationDomain("www-actar-get.cea.fr");
	    QApplication::setApplicationName(QObject::tr("cobo-frame-viewer"));
//#if QT_VERSION >= 0x040200
//	    app.setStyleSheet("QTabWidget#coboTabWidget { border-image: url(:/get/logo.png); }");
//	    app.setStyleSheet("QTabWidget { border-image: url(:/get/logo.png) }");
//#endif

	    LOG_INFO() << "Loading frame formats...";
		if (vm.count("format"))
		{
			try
			{
				LOG_DEBUG() << "Loading formats from file '" << fmtfile << '\'';
				mfm::FrameDictionary::instance().addFormats(fmtfile);
			}
			catch (const CCfg::Exception & e)
			{
				LOG_FATAL() << e.what();
				status = EXIT_FAILURE;
				return status;
			}
		}
		else
		{
			std::string formatPath = "CoboFormats.xcfg";
			try
			{
				LOG_DEBUG() << "Loading formats from file '" << formatPath << '\'';
				mfm::FrameDictionary::instance().addFormats(formatPath);
			}
			catch (const CCfg::Exception & e)
			{
				LOG_WARN() << e.what();
				try
				{
					formatPath = QString("%1/../share/get-bench/format/CoboFormats.xcfg").arg(QCoreApplication::applicationDirPath()).toStdString();
					LOG_DEBUG() << "Loading formats from file '" << formatPath << '\'';
					mfm::FrameDictionary::instance().addFormats(formatPath);
				}
				catch (const CCfg::Exception & e)
				{
					LOG_FATAL() << e.what() << "\nCould not load frame formats!";
				}
			}
		}

		// Theme
#if QT_VERSION >= 0x040600
		if (not QIcon::hasThemeIcon("application-exit"))
		{
			// If there is no default working icon theme (Windows, Mac OS X)
			// then we provide fallback theme
			QIcon::setThemeName("oxygen");
		}
#endif

		get::CoBoFrameViewer viewer;
		viewer.setChannels(channels);
		viewer.selectCurrentDisplay(agetIdx);
		if (vm.count("input-file"))
		{
			if (not viewer.setFrame(filename, frameIndex, byteOffset))
			{
				status = EXIT_FAILURE;
				return status;
			}
			if (vm.count("play"))
			{
				viewer.startAutoscan();
			}
		}
		if (vm.count("3d"))
		{
			viewer.set3dDisplayMode();
		}
		viewer.show();
		viewer.raise();
		status = app.exec();
	}
	catch (std::exception & e)
	{
		LOG_FATAL() << e.what();
		status = EXIT_FAILURE;
	}
	return status;
}

