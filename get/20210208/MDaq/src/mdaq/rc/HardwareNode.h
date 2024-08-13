/**
 * @file HardwareNode.h
 * @date 8 oct. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors:
 *
 * This file is part of the MDaq software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
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
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#ifndef mdaq_rc_HardwareNode_h_INCLUDED
#define mdaq_rc_HardwareNode_h_INCLUDED

#include "mdaq/hw/Control.h"
#include "mdaq/hw/RemoteRegisterProcessor.h"
#include "mdaq/DaqControl.h"
#include "utl/net/SocketAddress.h"

#include ICE_HEADER(Config.h)
#ifdef ICEE_INT_VERSION
#include ICE_HEADER(Handle.h)
#else
#include <IceUtil/Handle.h>
#endif

namespace mdaq
{
/**
 * Namespace for Run Control.
 */
namespace rc
{

//_________________________________________________________________________________________________
/**
 * A structure with all information related to an hardware node.
 */
struct HardwareNode : boost::noncopyable
{
	HardwareNode (
			::utl::net::SocketAddress endpoint = ::utl::net::SocketAddress(),
			 const hw::NodePrx& proxy = hw::NodePrx(0),
			 bool isBatched = false
	);
	operator bool() const;

	void updateDevices();
	hw::RegisterProcessor::Map::const_iterator findRegProc(const std::string & deviceName) const;
	hw::RegisterProcessor::Map::iterator findRegProc(const std::string & deviceName);
	void selectDevice(const std::string & deviceName);

	hw::RemoteRegisterProcessor& regProc(const std::string& deviceName);
	hw::RemoteRegisterProcessor& currentRegProc();
	hw::DevicePrx& currentDevicePrx();
	hw::DevicePrx& findDevice(const std::string & deviceName);
	std::string currentDeviceName();

	void connectDaqCtrl();
	mdaq::DaqCtrlNodePrx& daqCtrl();
	void connectRunProcessor(const std::string & dataRouterAddress, const std::string & dataRouterType);
	void disconnectRunProcessor();

	::utl::net::SocketAddress endpoint;
	hw::NodePrx hwServerCtrlPrx;
	DaqCtrlNodePrx dataSenderCtrlPrx;
	hw::RegisterProcessor::Map regProcMap; ///< List of remote register processors in node.
	hw::RegisterProcessor::Map::iterator currentRegProc_;

	const bool isBatched;
};

//_________________________________________________________________________________________________
} /* namespace rc */

typedef std::auto_ptr<rc::HardwareNode> HardwareNodePtr;

} /* namespace mdaq */
#endif /* mdaq_rc_HardwareNode_h_INCLUDED */
