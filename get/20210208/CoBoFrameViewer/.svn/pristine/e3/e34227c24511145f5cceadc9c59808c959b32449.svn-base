/*
 * $Id$
 * @file main.cpp
 * @created 10 juil. 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: patrick.sizun@cea.fr
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of
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

#include "graw2root.h"
#include <mfm/FrameDictionary.h>
#include <utl/Logging.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

int main(int argc, char* argv[])
{
	int status = EXIT_SUCCESS;
	std::string input_filename;
	std::string output_filename;
	std::string fmtfile;
	size_t frameOffset, frameCount;
	std::string config_filename;
	try
	{
		// Declare the supported options.
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message")
			("t2k", "produce ROOT tree of T2K (uncompressed) events")
			("input-file,i",   po::value<std::string>(&input_filename), "input raw frame filename")
			("output-file,o",   po::value<std::string>(&output_filename), "output ROOT filename")
			("format,f", po::value<std::string>(&fmtfile), "overrides format description file")
			("start,s", po::value<size_t>(&frameOffset)->default_value(0u), "index of first frame")
			("count,n", po::value<size_t>(&frameCount)->default_value(0u), "number of frames to convert (0 for all)")
			("config,c", po::value<std::string>(&config_filename), "path of configuration to embed");
		;
		po::positional_options_description p;
		p.add("input-file", 1);
		p.add("output-file", 1);

		// Parse options
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			LOG_INFO() << "Usage: " << argv[0] << " [--t2k] input.graw [output.root] [OPTIONS...]";
			LOG_INFO() << desc;
			return EXIT_FAILURE;
		}

		if (not vm.count("input-file"))
		{
			LOG_FATAL() << "No input file specified.\nUsage: " << argv[0] << " [--t2k] input.graw [output.root] [OPTIONS...]";
			return EXIT_FAILURE;
		}
		if (not vm.count("output-file"))
		{
#if BOOST_FILESYSTEM_VERSION >= 3
			std::string extension  = fs::path(input_filename).extension().string();
			if (".graw" == extension)
			{
				output_filename = fs::path(input_filename).filename().replace_extension(".root").string();
			}
			else
			{
				output_filename = fs::path(input_filename).filename().string() + std::string(".root");
			}
#else
			std::string extension  = fs::path(input_filename).extension();
			if (".graw" == extension)
			{
				output_filename = fs::path(input_filename).replace_extension(".root").filename();
			}
			else
			{
				output_filename = fs::path(input_filename).filename() + std::string(".root");
			}
#endif
		}

		// Load frame formats
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
				LOG_FATAL() << e.what() << "\nCould not load frame formats!";
				status = EXIT_FAILURE;
				return status;
			}
		}

		if (not vm.count("t2k"))
			status = get::graw2root(input_filename, output_filename, frameCount, frameOffset, config_filename);
		else
			status = get::graw2t2k(input_filename, output_filename, frameCount, frameOffset, config_filename);
	}
	catch(const std::exception & e)
	{
		LOG_FATAL() << e.what();
		status = EXIT_FAILURE;
	}
	return status;
}



