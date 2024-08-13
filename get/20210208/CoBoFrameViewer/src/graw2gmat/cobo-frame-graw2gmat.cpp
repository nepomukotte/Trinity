// $Id$
#include "get/CoBoEvent.h"
#include <CCfg/Exception.h>
#include <mfm/FrameDictionary.h>
#include <mfm/Frame.h>
#include <mfm/BitField.h>
#include <utl/Logging.h>
using namespace mfm;
using namespace get;

#include <QApplication>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
//#include <boost/timer.hpp>

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <vector>

bool write_matrix(const get::CoBoEvent & event, std::ofstream & out)
{
	size_t const maxAsics = 4u;
	size_t const maxChannels = 68u;
	size_t const maxSamples = 512u;
	try
	{
		for (size_t agetIdx=0; agetIdx < maxAsics; ++agetIdx)
		{
			for (size_t chanIdx = 0; chanIdx < maxChannels; chanIdx++)
			{
				const get::Channel & channel = event.channel(chanIdx, agetIdx);
				const uint16_t* samples = channel.data();
				out.write((const char*) samples, maxSamples * sizeof(uint16_t));
			}
		}
	}
	catch (const get::CoBoEvent::ChannelNotFound & e)
	{
		LOG_WARN() << e.what();
		LOG_ERROR() << "Error writing data to binary matrix file (missing channels)";
		return false;
	}
	return true;
}

bool convert_raw_frame(mfm::Frame & frame, std::ofstream & out)
{
	int status = EXIT_SUCCESS;

	try
	{
		get::CoBoEvent coboEvent;
		//coboEvent.fromFrame(frame);
		coboEvent.fromFrame_hardcoded(frame);
		bool ok = write_matrix(coboEvent, out);
		return ok;
	}
	catch (const mfm::Exception &e)
	{
		LOG_FATAL() << e.what();
		status = EXIT_FAILURE;
	}
	return (EXIT_SUCCESS == status);
}

bool graw2gmat(const std::string input_filename, const std::string & output_filename, const size_t maxFrames=0, const size_t frameOffset=0)
{
	LOG_DEBUG() << "Opening file '" << input_filename << "'...";
	std::ifstream in;
	in.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try
	{
		in.open(input_filename.c_str(), std::ios::in | std::ios::binary);
	}
	catch (const std::ifstream::failure & e)
	{
		LOG_ERROR() << "Could not open file '" << input_filename << "': " << strerror(errno);
		return false;
	}

	// Create output file
	LOG_INFO() << "Creating file '" << output_filename << "'...";
	std::ofstream out(output_filename.c_str(), std::ios::binary | std::ios::trunc);

	// Loop over frames in input file
	std::auto_ptr<Frame> frame;
	Frame::seekFrame(in, frameOffset);
	size_t frameCount = 0;
	bool ok=true;
	try
	{
		do
		{
			if (maxFrames > 0 and frameCount >= maxFrames) break;
			// Read frame
			try
			{
				//boost::timer timer;
				LOG_INFO() << "Converting frame #" << (frameOffset + frameCount++) << "...";
				frame = Frame::read(in);
				ok = convert_raw_frame(*frame.get(), out);
				//LOG_INFO() << "Frame converted in " << timer.elapsed() << " s";
				if (not ok) break;
			}
			catch (const std::ifstream::failure & e)
			{
				if (in.rdstate() & std::ifstream::eofbit)
				{
					LOG_WARN() << "EOF reached.";
					break;
				}
				else
				{
					LOG_ERROR() << "Error reading frame: " << e.what();
				}
				return false;
			}
		}
		while (in.good());
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << "Error: " << e.what();
		ok = false;
	}

	// Close input file
	in.close();

	// Close output file
	LOG_DEBUG() << "Closing output file";
	out.close();
	if (not ok)
	{
		LOG_INFO() << "Deleting output file";
		fs::remove(fs::path(output_filename));
	}

	return ok;
}

int main(int argc, char* argv[])
{
	int status = EXIT_SUCCESS;

	try
	{
		// List of options
		std::string input_filename;
		std::string output_filename;
		std::string fmtfile;
		size_t frameOffset, frameCount;

		// Define options
		po::options_description allowed_options("Options");
		allowed_options.add_options()("help,h", "help message")
    	     ("input-file,i",   po::value<std::string>(&input_filename), "input raw frame filename")
    	     ("output-file,o",   po::value<std::string>(&output_filename), "output matrix filename")
			 ("format,f", po::value<std::string>(&fmtfile), "overrides format description file")
			 ("start,s", po::value<size_t>(&frameOffset)->default_value(0u), "index of first frame")
			 ("count,n", po::value<size_t>(&frameCount)->default_value(0u), "number of frames to convert (0 for all)");
		po::positional_options_description p;
		p.add("input-file", 1);
		p.add("output-file", 1);

		// Parse options
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).
		          options(allowed_options).positional(p).run(), vm);
		po::notify(vm);
		// Process options
		if (vm.count("help") or not vm.count("input-file"))
		{
			LOG_INFO() << "Usage: " << argv[0] << " -i frames.graw [-o frames.gmat] [OPTIONS...]";
			LOG_INFO() << allowed_options << '\n';
			status = EXIT_FAILURE;
			return status;
		}
		if (not vm.count("output-file"))
		{
#if BOOST_FILESYSTEM_VERSION >= 3
			output_filename = fs::path(input_filename).filename().replace_extension(".gmat").string();
#else
			output_filename = fs::path(input_filename).replace_extension(".gmat").filename();
#endif
		}

		QApplication app(argc, argv);

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
					status = EXIT_FAILURE;
					return status;
				}
			}
		}

		LOG_INFO() << "Converting '" << input_filename << "' into '" << output_filename << "'...";
		graw2gmat(input_filename, output_filename, frameCount, frameOffset);
	}
	catch (std::exception & e)
	{
		LOG_FATAL() << e.what();
		status = EXIT_FAILURE;
	}
	return status;
}

