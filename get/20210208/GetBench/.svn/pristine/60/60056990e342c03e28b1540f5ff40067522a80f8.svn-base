/**
 * @file GetHwServer.h
 * Server that runs on the target device.
 * It permits to communicate with the target device through the MDaq interface
 */

#ifndef GetHwServer_INCLUDED
#define GetHwServer_INCLUDED

#include "get/hw/server/HardwareNodeImpl.h"
#include "mdaq/DefaultPortNums.h"
#include ICE_HEADER(ObjectAdapter.h)
#include ICE_HEADER(Communicator.h)
#include ICE_HEADER(Initialize.h)
#include ICE_HEADER(Properties.h)
#include <string>

/**
 * @brief Running and stopping the server on the PPC440 target device
 *
 * The server should be launched before loading the MDaq structure into the device
 */
class GetHwServer
{
public:
	static int  run(int argc, char* argv[]);
	static int  runInThread(int argc, char* argv[]);
	static void stop();
private:
	GetHwServer() {}
	static Ice::CommunicatorPtr ic;
};

#endif //GetHwServer_INCLUDED
