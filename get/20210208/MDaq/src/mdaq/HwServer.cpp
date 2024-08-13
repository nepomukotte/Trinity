/* =====================================================================================================================
 * HwServer.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * Helper functions for launching a Hardwarenode Server.
 * Created on: 12 February 2010 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
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

#include "mdaq/hw/server/HwNodeImpl.h"
#include "mdaq/DefaultPortNums.h"
#include "mdaq/utl/Server.h"
#include "utl/Logging.h"

#ifdef VXWORKS
#include <vxWorks.h>
#include <taskLib.h>
#endif

using ::mdaq::utl::Server;
using ::mdaq::hw::server::HwNodeImpl;
using ::mdaq::utl::CmdLineArgs;

static Server* hwServer = 0;

int hwServerRunOnAddressWithArgs(const utl::net::SocketAddress & sockAddr, const CmdLineArgs& args = CmdLineArgs())
{
	try
	{
		hwServer = &Server::create(sockAddr, args).addServant("HwNode", new HwNodeImpl).start();
#ifndef VXWORKS
		hwServer->ic()->waitForShutdown();
#endif
		return 0;
	}
	catch (const std::exception& e) { LOG_ERROR() << e.what() ; }
	catch (const std::string&  msg) { LOG_ERROR() <<   msg    ; }
	catch (const char*         msg) { LOG_ERROR() <<   msg    ; }
	catch (...)                     { LOG_ERROR() << "UNKNOWN"; }
	return 1;
}

int hwServerRunOnPortWithArgs(int portNum, const CmdLineArgs& args = CmdLineArgs())
{
	return hwServerRunOnAddressWithArgs(utl::net::SocketAddress(0, portNum), args);
}

int hwServerRunOnPort(int portNum = ::mdaq::Default::hwServerCtrlPortNum)
{
	return hwServerRunOnPortWithArgs(portNum);
}

int hwServerRun()
{
	return hwServerRunOnPort();
}

int hwServerStop()
{
	try
	{
		if (hwServer != 0)
		{
			hwServer->stop();
		}
		return 0;
	}
	catch (const std::exception& e) { LOG_ERROR() << e.what() ; }
	catch (const std::string&  msg) { LOG_ERROR() <<   msg    ; }
	catch (const char* msg)         { LOG_ERROR() <<   msg    ; }
	catch (...)                     { LOG_ERROR() << "UNKNOWN"; }
	return 1;
}

#ifdef VXWORKS

namespace mdaq
{
class Launcher
{
	Launcher()
	{

		// The server is launched using taskspawn() to prevent it from inheriting priority 1 from shell.
		if (taskSpawn("HwServer", 10, VX_FP_TASK, 0x10000, (FUNCPTR) hwServerRun, 0,1,2,3,4,5,6,7,8,9) == ERROR)
		{
			LOG_ERROR() << "Failed to launch Hardware Server";
		}
	}

	static Launcher launcher;
};

Launcher Launcher::launcher;
}

#else

int main(int argc, char* argv[])
{
	CmdLineArgs args = CmdLineArgs::build(argc, argv);

	// Set ICE communicator initialization data
	::mdaq::utl::Server::setInitData(args);

	// Check if address was specified
	if (args.size() > 1)
	{
		::utl::net::SocketAddress sockAddr(args[1]);
		if (sockAddr.isValid())
		{
			if (sockAddr.port() == 0)
			{
				sockAddr.port() = ::mdaq::Default::hwServerCtrlPortNum;
			}
			return hwServerRunOnAddressWithArgs(sockAddr, args);
		}
		else
		{
			LOG_ERROR() << "Invalid IP address / port: " << sockAddr;
			return EXIT_FAILURE;
		}
	}

	return hwServerRunOnPortWithArgs(::mdaq::Default::hwServerCtrlPortNum, args);
}
#endif

