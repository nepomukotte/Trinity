/*
 * hwLoad.cpp
 *
 *  Created on: Mar 17, 2010
 *      Author: jchavas
 */

/**
 * @file hwLoad.cpp
 * It loads into a target device a xml CompoundConfig configuration file
 */

#include "mdaq/cfg/NodeAtTop.h"
#include "mdaq/DefaultPortNums.h"
#include "CCfg/Io/Document.h"
#include <iostream>
#include <cstdlib>
#include <unistd.h>

int main(int argc, char **argv)
{
	/****************************************
	 * Intepretation of the command lines
	 ***************************************/

	if ( (argc > 1) && (strcmp(argv[1],"-h") == 0) )
	{
		std::cout << "Usage: " << argv[0] << " ip file" << std::endl;
		std::cout << "Usage: " << argv[0] << " ip:port file" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58 CCfg_file.xml" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58:10001 CCfg_file.xml" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58:10001 CCfg_file.xml --Ice.Trace.Network=3 --Ice.Trace.Protocol=3" << std::endl;
		return 0;
	}

	if (argc < 3)
	{
		std::cout << "*** there must be 2 parameters, the IP of the target and the name of the configuration file"  << std::endl;
		std::cout << "Usage: " << argv[0] << " ip file" << std::endl;
		std::cout << "Usage: " << argv[0] << " ip:port file" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58 CCfg_file.xml" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58:10001 CCfg_file.xml" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58:10001 CCfg_file.xml --Ice.Trace.Network=3 --Ice.Trace.Protocol=3" << std::endl;
		return 1;
	}

	std::string sfile_cfg(argv[2]);

	/****************************************
	 * One describes a structure NoteAtTop : fullCfg.
	 * Nothing is loaded to the remote target so far
	 ***************************************/

	mdaq::cfg::NodeAtTop fullCfg;

	try
	{
		CCfg::Io::Document doc(sfile_cfg);
		CCfg::CConfig cfg = doc.getConfig();
		fullCfg.fromCCfg(cfg);
	}
	catch (const CCfg::Io::Document::IoError& err)
	{
		std::cout << err.what() << std::endl;
		return 1;
	}
	catch (const std::exception& err)
	{
		std::cout << "*** " << err.what() << std::endl;
		return 1;
	}
	catch (const char *msg)
	{
		std::cout << msg << std::endl;
		return 1;
	}

	/****************************************
	 * Ice stuff
	 ***************************************/

	Ice::CommunicatorPtr ic;
	mdaq::hw::NodePrx node;
	std::string strProxy("HwNode:default");
	{
		std::string str_argv1(argv[1]);
		size_t found;

		strProxy += " -h ";
		found = str_argv1.find(':',0);
		strProxy += (found != std::string::npos) ? str_argv1.substr(0,found) : str_argv1;

		strProxy += " -p ";
		strProxy += (found != std::string::npos) ? str_argv1.substr(found+1,str_argv1.length()+1) : mdaq::Default::hwServerCtrlPortStr;
	}

	try
	{
		ic = Ice::initialize(argc, argv);
		Ice::ObjectPrx base = ic->stringToProxy(strProxy);
		node = mdaq::hw::NodePrx::checkedCast(base);
		if(not node) throw "Invalid proxy";
	}
	catch (const Ice::Exception& ex)
	{
		std::cerr << "*** " << ex.what() << std::endl;
		return 1;
	}
	catch (const char* msg)
	{
		std::cerr << msg << std::endl;
		return 1;
	}

	/****************************************************
	 * loads the description to the target and instantiates the node devices
	 * It first destroys the configuration loaded on the target
	 ****************************************************/
	mdaq::hw::DevicePrx devicePrx;

	std::cout << "Reloading NOW the configuration file " << sfile_cfg << std::endl;
	node->destroy();
	node->create(fullCfg);

	// takes the index of the first expected device
	std::string strDevice = fullCfg.devices[0].descr.name;
	try
	{
		devicePrx = node->findDevice(strDevice);
	}
	catch (const ::mdaq::utl::CmdException& e)
	{
		std::cout << "ERROR : the configuration file has not been loaded properly" << std::endl;
		std::cout << "It has been a HardwareNotFound exception" << std::endl;
	}


	return 0;
}
