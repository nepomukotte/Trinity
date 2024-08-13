/*
 * $Id: SocketAddress.h 1338 2013-10-17 07:53:32Z psizun $
 * @file SocketAddress.h
 * @created 11 mai 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar, Patrick Sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the Utilities package of the Mordicus Real-Time
 * Software Development Framework.
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

#ifndef utl_net_SocketAddress_h_INCLUDED
#define utl_net_SocketAddress_h_INCLUDED

#include <utl/net/IpAddress.h>
#include <string>

namespace utl {
namespace net {
//__________________________________________________________________________________________________
/**
 * An IP socket address, defined as a combination of an IP Internet address and a port number.
 */
class SocketAddress
{
public:
	SocketAddress(uint16_t port = 0);
	SocketAddress(const utl::net::IpAddress &, uint16_t port);
	SocketAddress(in_addr_t address, uint16_t port);
	SocketAddress(const std::string& dotDecimalIpString, uint16_t port);
	SocketAddress(const std::string& dotDecimalIpAndPortString);
	~SocketAddress() {}

	const utl::net::IpAddress & ipAddress() const { return ipAddress_; }
	      utl::net::IpAddress & ipAddress()       { return ipAddress_; }
	uint16_t   port() const { return port_; }
	uint16_t & port()       { return port_; }

	bool isValid() const { return ipAddress().isValid(); }

	std::string iceEndPointString(const std::string & protocol = "default") const;
	std::string iceProxyString(const std::string & identity, const std::string & protocol = "default", const std::string options="") const;
	std::string toString() const;
	void setPort(const std::string & port);
	SocketAddress& fromString(const std::string& decimalString);
	SocketAddress& from_sockaddr_in(const struct sockaddr_in& bsd_sockaddr);
	friend std::ostream& operator<<(std::ostream& stream, const SocketAddress& url);

private:
	utl::net::IpAddress ipAddress_; ///< Internet IPv4 address
	uint16_t port_; ///< 16 bit port number
};
//__________________________________________________________________________________________________
} /* namespace net */
} /* namespace utl */
#endif /* utl_net_SocketAddress_h_INCLUDED */
