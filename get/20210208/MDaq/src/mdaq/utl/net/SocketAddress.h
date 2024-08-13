/*
 * $Id: SocketAddress.h 834 2012-09-11 16:20:02Z anvar $
 * @file SocketAddress.h
 * @created 11 May 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::utl::net::SocketAddress
 * Created on: 11 May, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Patrick Sizun (patrick.sizun@cea.fr)
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

#ifndef mdaq_utl_net_SocketAddress_h_INCLUDED
#define mdaq_utl_net_SocketAddress_h_INCLUDED

#include <utl/net/SocketAddress.h>

// Forward declarations
namespace CCfg
{
	class CConfig;
}

namespace mdaq {
namespace utl {
namespace net {
//__________________________________________________________________________________________________
/**
 * A combination of an IP address and a port number, which can be saved or loaded to or from a configuration.
 */
class SocketAddress : public ::utl::net::SocketAddress
{
public:
	SocketAddress(uint16_t port = 0);
	SocketAddress(in_addr_t address, uint16_t port);
	SocketAddress(const std::string& dotDecimalIpString, uint16_t port);
	SocketAddress(const std::string& dotDecimalIpAndPortString);
	SocketAddress(CCfg::CConfig& attr);
	SocketAddress & fromCConfig(const CCfg::CConfig& obj);
	void toCConfig(const CCfg::CConfig& obj) const;
};
//__________________________________________________________________________________________________
} /* namespace net */
} /* namespace utl */
} /* namespace mdaq */
#endif /* mdaq_utl_net_SocketAddress_h_INCLUDED */
