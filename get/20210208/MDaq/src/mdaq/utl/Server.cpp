/* =====================================================================================================================
 * Server.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::utl::Server
 * Created on: Jan 5, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the MDaq data acquisition framework @ CEA-Irfu.
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

#include "Server.h"

#ifdef VXWORKS
#include <vxWorks.h>
#include <hostLib.h>
#endif

#include "utl/Logging.h"

using ::utl::net::SocketAddress;

namespace mdaq
{
namespace utl
{

::Ice::CommunicatorPtr Server::ic_; ///< The ICE communicator runtime.
CmdLineArgs Server::initData_; ///< Initialization data for the ICE communicator.
Server::Map Server::servers_; ///< The servers for every endpoint.

/**
 * Returns the instance of the Server listening to specified port, if it has been created. If server endpoint is
 * different from the one already created, the server is re-created with the new endpoint.
 * @return If created, the reference to the Server instance.
 * @throw NotCreatedException if Server has not been created (instantiated).
 */
Server& Server::create(const SocketAddress& ip, const CmdLineArgs& args)
{
	SocketAddress normIp = normalized(ip);
	Map::iterator it = servers_.find(normIp.toString());
	if (it != servers_.end())
	{
		delete servers_[normIp.toString()];
	}
	servers_[normIp.toString()] = new Server(normIp, args);
	return *servers_[normIp.toString()];
}

/**
 * Compute normalized IP address & port specification from a given one. Namely, if the address field of the given ip is
 * zero, it computes a default address (which is different on various platforms).
 * @param ip The given IP address & port specification.
 * @return The normalized IP address & port specification.
 */
SocketAddress Server::normalized(const SocketAddress& ip)
{
	SocketAddress localIp(ip);
#ifdef VXWORKS
	if (ip.ipAddress().in_addr() == 0) // No host defined
	{

		// Fabricating local endpoint for the VxWorks case. This feature needs a tailored VxWorks image in which the
		// default hostname vxTarget is properly associated with the target's IP interface.
		localIp.ipAddress().in_addr() = hostGetByName("vxTarget");
	}
#endif
	return localIp;
}

/**
 * Finds server listening to (normalized) specified ip.
 * @param ip The IP address & port specifications of the server.
 * @return A reference to the server if found.
 * @throw Server::NotCreatedException if server not found.
 */
Server& Server::find(const SocketAddress& ip)
{
	if (servers_.find(normalized(ip).toString()) == servers_.end())
		throw NotCreatedException(normalized(ip));

	return *servers_[normalized(ip).toString()];
}

/**
 * Destroys server with specified address.
 * @param addr The IP address & port specifications of the server.
 */
void Server::destroy(const ::utl::net::SocketAddress& addr)
{
	const std::string key = normalized(addr).toString();
	if (servers_.find(key) != servers_.end())
	{
		delete servers_[key];
		servers_.erase(key);
	}
}

/**
 * Server constructor. Needs a name and a port number.
 * @param ip The IP address & port (endpoint) on which the server is to listen for requests.
 */
Server::Server(const SocketAddress& ip, const CmdLineArgs& args) : isStarted_(false), ip_(ip)
{
	setInitData(args);
}

/**
 * Server destructor. Stops the server before returning from destructor.
 */
Server::~Server()
{
	stop();
}

/**
 * Sets the initialization data for the ICE communicator.
 * Does nothing of the ICE communicator has already been created.
 *
 * @param argc Number of arguments
 * @param argv Argument array.
 */
void Server::setInitData(int argc, char** argv)
{
	setInitData(CmdLineArgs::build(argc, argv));
}

/**
 * Sets the initialization data for the ICE communicator.
 * Does nothing of the ICE communicator has already been created.
 *
 * @param args Initialization data.
 */
void Server::setInitData(const CmdLineArgs & args)
{
	// Set initialization data
	if (not args.empty())
	{
		if (ic_ and initData_ != args)
		{
			LOG_WARN() << "ICE communicator has already been created, new initialization data will be ignored.";
		}
		else if (not ic_)
		{
			initData_ = args;
		}
	}
}

/**
 * Creates Ice communicator if not yet created.
 * @param newArgs Initialization data for communicator.
 * @return Smart pointer to the Ice communicator.
 */
Ice::CommunicatorPtr Server::ic(const CmdLineArgs & args)
{
	setInitData(args);

	// Create communicator
	if (!ic_)
	{
		// Ice::initialize removes Ice-specific command-line options from the argument vector.
		// We need to keep a copy of the initial argument vector for later.
		CmdLineArgs initData(initData_);
		ic_ = Ice::initialize(initData);
	}
	return ic_;
}

/**
 * Creates Ice object adapter if not yet created.
 * @return Smart pointer to the Ice object adapter.
 */
Ice::ObjectAdapterPtr Server::adapter()
{
	if (!adapter_)
	{
		const std::string endpoint = ip().iceEndPointString();
		LOG_DEBUG () << "Creating object adapter @ " << endpoint;
		adapter_ = ic()->createObjectAdapterWithEndpoints (std::string("Adapter") + ip().toString(), endpoint);
	}
	return adapter_;
}

/**
 * Starts the server by activating its object adapter. Automatically allocates communicator and adapter if not yet done.
 * Has no effect if the server is already started.
 * @return Reference to server object so as to allow chained calls.
 */
Server& Server::start()
{
	if (isStopped())
	{
		adapter()->activate();
		isStarted_ = true;
		MSG() << "STARTED server on " << ip();
	}
	else
	{
		LOG_INFO() << "Server ALREADY STARTED on " << ip();
	}
	return *this;
}

/**
 * Requests Server to stop (calls adapter deactivation but does not wait until done).
 * @return Reference to server object so as to allow chained calls.
 */
Server& Server::askForStop()
{
	if (isStarted())
	{
		adapter()->deactivate();
		isStarted_ = false;
	}
	else
	{
		LOG_DEBUG() << "Server ALREADY STOPPED on " << ip();
	}
	return *this;
}

/**
 * Waits for Server to stop (waits for adapter deactivation and then destroys adapter).
 * @return Reference to server object so as to allow chained calls.
 */
Server& Server::waitForStop()
{
	if (adapter_)
	{
		adapter()->waitForDeactivate();
		adapter()->destroy();
		adapter_ = 0;
		LOG_DEBUG() << "STOPPED server on " << ip();
	}
	return *this;
}

/**
 * @fn Server& Server::stop()
 * Stops the server (puts adapter on hold. This stops processing requests. If the server is already on hold, nothing
 * happens.
 * @return Reference to server object so as to allow chained calls.
 */

/**
 * Adds a new servant to this server.
 * @param servantName Servant identity within this server.
 * @param servant The servant object.
 * @return  Reference to server object so as to allow chained calls.
 * @throw Ice::AlreadyRegisteredException if a servant has already been registered under this name.
 */
Server& Server::addServant(const std::string& servantName, Ice::ObjectPtr servant)
{
	LOG_DEBUG() << "Adding servant with id: " << servantName;
	adapter()->add(servant, ic()->stringToIdentity(servantName));
	return *this;
}

/**
 * Adds a new servant to this server.
 * @param servant The servant object.
 * @return Identity name given to the servant.
 */
std::string Server::addServantWithUUID(Ice::ObjectPtr servant)
{
	LOG_DEBUG() << "Adding servant with UUID ...";
	Ice::ObjectPrx proxy = adapter()->addWithUUID(servant);
	return proxy->ice_getIdentity().name;
}

Ice::ObjectPrx Server::servantPrx(const std::string& servantName)
{
	return ic()->stringToProxy(ip_.iceProxyString(servantName));
}

} // namespace utl
} // namespace mdaq
