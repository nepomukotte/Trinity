/**
 * @file cobo-frame-merge.cpp
 * @date Oct 10, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: cobo-frame-merge.cpp 1508 2015-10-15 16:24:29Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the CoBoFrameViewer software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
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
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include <utl/Logging.h>
#include <mfm/FrameDictionary.h>
#include "CoBoFrameMerger.h"

#include <QCoreApplication>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/cstdint.hpp>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

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

std::string build_outputFileName(const std::string & inputFileName)
{
	std::string fileName;
#if BOOST_FILESYSTEM_VERSION >= 3
	std::string extension  = fs::path(inputFileName).extension().string();
	if (".graw" == extension)
	{
		fileName = fs::path(inputFileName).filename().replace_extension(".gmer").string();
	}
	else
	{
		fileName = fs::path(inputFileName).filename().string() + std::string(".gmer");
	}
#else
	std::string extension  = fs::path(inputFileName).extension();
	if (".graw" == extension)
	{
		fileName = fs::path(inputFileName).replace_extension(".gmer").filename();
	}
	else
	{
		fileName = fs::path(inputFileName).filename() + std::string(".gmer");
	}
#endif
	return fileName;
}

int main(int argc, char* argv[])
{
	int status = EXIT_SUCCESS;
	std::string fmtfile;
	std::string output_filename;
	size_t maxMergedEvents = 0;
	uint32_t minEventID = 0;
	uint32_t maxEventID = 0;
	uint32_t maxFileSize_MiB = 0;
	try
	{
		// Declare the supported options.
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message")
			("format,f", po::value<std::string>(&fmtfile), "overrides format description file")
			("input-file,i", po::value< std::vector< std::string > >(), "input raw MFM files")
			("output-file,o", po::value<std::string>(&output_filename), "name of output merged data file")
			("max-merged-events,n", po::value< size_t >(&maxMergedEvents)->default_value(0), "max number of merged events")
			("min-event-id,m", po::value< uint32_t >(&minEventID)->default_value(0), "minimum event identifier")
			("max-event-id,M", po::value< uint32_t >(&maxEventID)->default_value(UINT32_C(4294967295)), "maximum event identifier")
			("max-file-size,s", po::value< uint32_t >(&maxFileSize_MiB)->default_value(UINT32_C(1024)), "maximum size of each output file [MiB] or 0")
		;
		po::positional_options_description p;
		p.add("input-file", -1);

		// Parse options
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			LOG_INFO() << "Usage: " << argv[0] << " input.graw [input2.graw ...] [OPTIONS...]";
			LOG_INFO() << desc;
			return EXIT_FAILURE;
		}

		// Input files
		if (not vm.count("input-file"))
		{
			LOG_FATAL() << "No input file specified.\n";
			LOG_INFO() << desc;
			return EXIT_FAILURE;
		}
		std::vector < std::string > input_filenames = vm["input-file"].as< std::vector< std::string> >();

		// Output file
		if (not vm.count("output-file"))
		{
			if (input_filenames.size() >= 2)
			{
				LOG_FATAL() << "No output file specified.\n";
				LOG_INFO() << desc;
				return EXIT_FAILURE;
			}
			else
			{
				output_filename = build_outputFileName(input_filenames.at(0));
			}
		}

		// Frame format
		QCoreApplication app(argc, argv);
	    QCoreApplication::setOrganizationName("actar-get");
	    QCoreApplication::setOrganizationDomain("www-actar-get.cea.fr");
	    QCoreApplication::setApplicationName(QObject::tr("cobo-frame-merge"));
		if (EXIT_SUCCESS != (status = load_format(fmtfile)))
			return status;

		// Map frames
		get::data::CoBoFrameMerger merger;
		for (size_t i=0; i < input_filenames.size(); ++i)
		{
			merger.addFrames(input_filenames[i]);
		}
		merger.print();

		// Merge
		merger.mergeByEventId(output_filename, minEventID, maxEventID, maxMergedEvents, maxFileSize_MiB*UINT64_C(0x100000));

	}
	catch(const std::exception & e)
    {
		LOG_FATAL() << e.what();
		status = EXIT_FAILURE;
    }
	return status;
}



