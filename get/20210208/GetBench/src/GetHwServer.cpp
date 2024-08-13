/* =====================================================================================================================
 * GetHwServer.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * Helper functions for launching a hardware and acquisition Server.
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

#ifdef VXWORKS
#include "get/daq/DaqCtrlNodeI.h"
#include <vxWorks.h>
#include <hostLib.h>
#endif

#include "get/cobo/MutantLinkManagerI.h"
#include "get/cobo/CtrlNodeI.h"
#include "get/mt/AlarmServiceI.h"
#include "mdaq/hw/server/HwNodeImpl.h"
#include "mdaq/DefaultPortNums.h"
#include "mdaq/utl/Server.h"
#include "utl/net/SocketAddress.h"
#include <cstring> // strerror
#include <iostream>

#ifdef VXWORKS
using ::get::daq::DaqCtrlNodeI;
#endif

using ::mdaq::utl::Server;
using ::mdaq::hw::server::HwNodeImpl;
using ::mdaq::utl::CmdLineArgs;
using ::utl::net::SocketAddress;

static Server* getHwServer = 0;

#ifdef VXWORKS
static Server* daqHwServer = 0;
#endif

int getHwServerRunOnAddressAndPortWithArgs(const SocketAddress & endpoint, const CmdLineArgs& args = CmdLineArgs())
{
	try
	{
		get::cobo::CtrlNodeI* ctrlNodePtr = new get::cobo::CtrlNodeI;

		getHwServer = &Server::create(endpoint, args)
			.addServant("HwNode", ctrlNodePtr)
			.addServant("MutantLinkManager", new get::cobo::MutantLinkManagerI)
			.addServant("AlarmService", new get::mt::AlarmServiceI(*ctrlNodePtr))
			.start();

#ifdef VXWORKS
		SocketAddress daqCtrlEndpoint(endpoint.ipAddress(), ::mdaq::Default::dataSenderCtrlPortNum);
		daqHwServer = &Server::create(daqCtrlEndpoint, args)
			.addServant("DaqCtrlNode" , new DaqCtrlNodeI(*ctrlNodePtr))
			.start();
#else
		getHwServer->waitForStop();
#endif
		return 0;
	}
	catch (const Ice::SocketException & e)
	{
		std::cerr << "*** Socket exception: " <<  std::strerror(e.error) << '.' << std::endl;
	}
	catch (const std::exception& e) { std::cerr << "*** " << e.what()  << std::endl; }
	catch (const std::string&  msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (const char*         msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (...)                     { std::cerr << "*** " << "UNKNOWN" << std::endl; }
	return 1;
}

int getHwServerRunOnPortWithArgs(int portNum, const CmdLineArgs& args = CmdLineArgs())
{
	return getHwServerRunOnAddressAndPortWithArgs(SocketAddress(INADDR_ANY, portNum), args);
}

int getHwServerRunOnPort(int portNum = ::mdaq::Default::hwServerCtrlPortNum)
{
	return getHwServerRunOnPortWithArgs(portNum);
}

int getHwServerRun()
{
	SocketAddress endpoint(INADDR_ANY, ::mdaq::Default::hwServerCtrlPortNum);

	// Get server IP address from target host name
	// R-CoBo: 1 interface mdaqNode alias vxTarget
	// M-CoBo: 2 interfaces CoBoCtrlNode and CoBoDataNode (1 has alias vxTarget)
#ifdef VXWORKS
	if (ERROR == (endpoint.ipAddress().in_addr() = hostGetByName("CoBoCtrlNode")))
	{
		endpoint.ipAddress().in_addr() = hostGetByName("vxTarget");
	}
#endif

	return getHwServerRunOnAddressAndPortWithArgs(endpoint);
}

int getHwServerStop()
{
	try
	{
		getHwServer->stop();
#ifdef VXWORKS
		daqHwServer->stop();
#endif
		return 0;
	}
	catch (const std::exception& e) { std::cerr << "*** " << e.what()  << std::endl; }
	catch (const std::string&  msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (const char*         msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (...)                     { std::cerr << "*** " << "UNKNOWN" << std::endl; }
	return 1;
}

#ifdef VXWORKS

namespace get
{
class Launcher
{
	Launcher()
	{
		if (taskSpawn("GetHwServer", 10, VX_FP_TASK, 0x10000, (FUNCPTR) getHwServerRun, 0,1,2,3,4,5,6,7,8,9) == ERROR)
		{
			fprintf(stderr, "*** Failed to launch GET Hardware Server1\n");
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
	return getHwServerRunOnPortWithArgs(::mdaq::Default::hwServerCtrlPortNum, args);
}
#endif

