/*
 * FullConfig.cpp
 *
 *  Created on: Mar 17, 2010
 *      Author: jchavas
 */

/**
 * @file hwFieldRead.cpp
 * It reads a value from a bitfield from a target device
 */

#include "mdaq/cfg/NodeAtTop.h"
#include "mdaq/DefaultPortNums.h"
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
		std::cout << "Usage: " << argv[0] << " ip device reg bitfield" << std::endl;
		std::cout << "Usage: " << argv[0] << " ip:port device reg bitfield" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58 AsAd sw_status s0" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58:10001 AsAd sw_status s0" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58:10001 AsAd sw_status s0 --Ice.Trace.Network=3 --Ice.Trace.Protocol=3" << std::endl;
		return 0;
	}

	if (argc < 5)
	{
		std::cout << "ERR  : there must be 4 parameters"  << std::endl;
		std::cout << "Usage: " << argv[0] << " ip device reg bitfield" << std::endl;
		std::cout << "Usage: " << argv[0] << " ip:port device reg bitfield" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58 AsAd sw_status s0" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58:10001 AsAd sw_status s0" << std::endl;
		std::cout << "Usage: " << argv[0] << " 132.166.37.58:10001 AsAd sw_status s0 --Ice.Trace.Network=3 --Ice.Trace.Protocol=3" << std::endl;
		return -1;
	}

	/****************************************
	 * Ice stuff
	 ***************************************/

	Ice::CommunicatorPtr ic;
	mdaq::hw::NodePrx node;
	std::string strProxy("HwNode:default");
	{
		strProxy += " -h ";
		std::string str_argv1(argv[1]);
		size_t found;
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
		std::cerr << "ERR: " << ex.what() << std::endl;
		return -1;
	}
	catch (const char* msg)
	{
		std::cerr << msg << std::endl;
		return -1;
	}

	/****************************************************
	 * checks if the Device is present on the remote target
	 ****************************************************/
	mdaq::hw::DevicePrx devicePrx;
	std::string strDevice = argv[2];

	try
	{
		devicePrx = node->findDevice(strDevice);
	}
	catch (const ::mdaq::utl::CmdException& e)
	{
		std::cerr << "ERR: Device[" << strDevice << "] is NOT present on the remote target" << std::endl;
		return -1;
	}

	/****************************************************
	 * user program
	 ****************************************************/

	int valField;
	try
	{
		valField = devicePrx->readField(argv[3], argv[4]);
		std::cout << valField << std::endl;
	}
	catch (const ::mdaq::utl::CmdException& e)
	{
		std::cerr << "ERR: Device[" << strDevice << "].Register[" << argv[3] << "].BitField["<< argv[4] << "] is NOT present on the remote target" << std::endl;
		return -1;
	}

	return valField;
}
