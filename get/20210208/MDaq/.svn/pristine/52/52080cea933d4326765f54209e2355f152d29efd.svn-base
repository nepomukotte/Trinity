/* =====================================================================================================================
 * Server.h
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

#ifndef mdaq_utl_Server_h_INCLUDED
#define mdaq_utl_Server_h_INCLUDED

#include "utl/net/SocketAddress.h"
#include "CmdLineArgs.h"
#include "utl/Exception.h"
#include ICE_HEADER(ICE_EDITION.h)
#include ICE_HEADER(Communicator.h)
#include ICE_HEADER(Object.h)
#include <boost/cstdint.hpp>
#include <map>

namespace mdaq
{
namespace utl
{

class Server
{
public:
	/**
	 * Exception thrown when Server is accessed before being created.
	 */
	class NotCreatedException : public ::utl::Exception
	{
	public:
		NotCreatedException(const ::utl::net::SocketAddress& ip) : index(Server::normalized(ip).toString()) {}
		NotCreatedException(uint16_t portNum) : index(Server::normalized(::utl::net::SocketAddress(0,portNum)).toString()) {}
		virtual ~NotCreatedException() throw() {}
		void writeMessage(std::ostream& stream) const { stream << "Server " << index << " not created"; }
		const std::string index;
	};

	/**
	 * This class extends a map by adding cleanup of the contained heap-allocated pointers.
	 */
	struct Map : public std::map<std::string, Server*>
	{
		virtual ~Map() { for (iterator it=begin(); it != end(); ++it) delete it->second; }
	};

public:
	static Server& create(const ::utl::net::SocketAddress& ip, const CmdLineArgs& args = CmdLineArgs());
	static Server& create(uint16_t portNum, const CmdLineArgs& args = CmdLineArgs()) { return create(::utl::net::SocketAddress(INADDR_ANY,portNum), args); }
	static ::utl::net::SocketAddress normalized(const ::utl::net::SocketAddress& ip);
	static Server& find(const ::utl::net::SocketAddress& ip);
	static Server& find(uint16_t portNum) { return find(::utl::net::SocketAddress(0,portNum)); }
	static void destroy(const ::utl::net::SocketAddress& addr);
	static void destroy(uint16_t portNum) { return destroy(::utl::net::SocketAddress(0,portNum)); }
	static Ice::CommunicatorPtr ic(const CmdLineArgs& args = CmdLineArgs());
	static void setInitData(int argc, char** argv);
	static void setInitData(const CmdLineArgs & newArgs = CmdLineArgs());
	static CmdLineArgs& initData() { return initData_; }

protected:
	Server(const ::utl::net::SocketAddress& ip, const CmdLineArgs& args);
	virtual ~Server();

public:
	Ice::ObjectAdapterPtr adapter();
	Server& start();
	Server& askForStop();
	Server& waitForStop();
	Server& stop() { return askForStop().waitForStop(); }
	Server& addServant(const std::string& servantName, Ice::ObjectPtr servant);
	std::string addServantWithUUID(Ice::ObjectPtr servant);
	Ice::ObjectPrx servantPrx(const std::string& servantName);
	::utl::net::SocketAddress ip() { return ip_; }
	bool isStarted() { return isStarted_; }
	bool isStopped() { return not isStarted_; }

private:
	static Ice::CommunicatorPtr ic_; ///< The ICE communicator runtime.
	static CmdLineArgs initData_; ///< Initialization data for the ICE communicator.
	static Map servers_; ///< The servers for every endpoint.

private:
	bool isStarted_;
	Ice::ObjectAdapterPtr adapter_;
	::utl::net::SocketAddress ip_; ///< The ip address and port number on which the server is listening (endpoint).
};


} // namespace utl
} // namespace mdaq

#endif /* mdaq_utl_Server_h_INCLUDED */
