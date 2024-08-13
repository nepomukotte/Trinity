/**
 * @file IpDetect.cpp
 * @created Nov 26, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::utl::net::IpDetect
 * Created on: Nov 26, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of Mordicus Real-Time Software Development Framework (MDaq).
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

#include "IpDetect.h"
#include <iostream>

#ifdef VXWORKS

#include <vxWorks.h>
#include <hostLib.h>
#include <netLib.h>
#include <inetLib.h>

#else

#include <unistd.h>
#include <net/if.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/ioctl.h>
#include <string.h>
#include <errno.h>

#endif

namespace mdaq
{
namespace utl
{
namespace net
{

/**
 * Obtains, in  network order number (in_addr_t), all IP addresses attached to
 * host on which it is instantiated.
 */
IpDetect::IpDetect()
{
#ifdef VXWORKS

	char hostName[128];
	gethostname(hostName, 128);

	//std::cout << "hostname = " << hostName << std::endl;
	addrinfo* ai;
	getaddrinfo(hostName, 0, 0, &ai);

	for (int i=0; ai != 0; ++i)
	{
		push_back((((sockaddr_in*)ai->ai_addr)->sin_addr).s_addr);
		ai = ai->ai_next;
	}

	freeaddrinfo(ai);

#else

	int iSocket = -1;
	// create a socket to be used when calling ioctl().
	if ((iSocket = socket(PF_INET, SOCK_DGRAM, 0)) < 0)
	{
		throw std::string(strerror(errno));
	}

	// if_nameindex - return all network interface names and indexes
	struct if_nameindex* pIndex = if_nameindex();

	while ((pIndex != NULL) and (pIndex->if_name != NULL))
	{
		struct ifreq req;

		strncpy(req.ifr_name, pIndex->if_name, IFNAMSIZ);

		// ioctl - control a STREAMS device
		if (ioctl(iSocket, SIOCGIFADDR, &req) < 0)
		{
			if (errno == EADDRNOTAVAIL)
			{
				// This interface has no assigned IP address: try the other ones.
			}
			else
			{
				close(iSocket);
				throw std::string(strerror(errno));
			}
		}
		else
		{
			push_back((((struct sockaddr_in*) &req.ifr_addr)->sin_addr).s_addr);
		}
		++pIndex;
	}

#endif

}

/**
 * Checks whether specified IP address belongs to this host.
 * @param ipAddrNum in binary network order (big-endian) format.
 * @return true if parameter corresponds to one of the IP addresses detected on
 * the NIC of this host.
 */
bool IpDetect::belongs(in_addr_t ipAddrNum) const
{
	for (const_iterator it = begin(); it != end(); ++it) {
		if (*it == ipAddrNum) {
			return true;
		}
	}
	return false;
}

/**
 * Checks whether specified IP address belongs to this host.
 * @param ipDotDec in dot-decimal notation.
 * @return true if parameter corresponds to one of the IP addresses detected on
 * the NIC of this host.
 */
bool IpDetect::belongs(const std::string& ipDotDec) const
{
	in_addr ipNet;
#ifdef VXWORKS
	inet_aton(const_cast<char*>(ipDotDec.c_str()), &ipNet);
#else
	inet_aton(ipDotDec.c_str(), &ipNet);
#endif
	return belongs(ipNet.s_addr);
}

} // namespace net
} // namespace utl
} // namespace mdaq

