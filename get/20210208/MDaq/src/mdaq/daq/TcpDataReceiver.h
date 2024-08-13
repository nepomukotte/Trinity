/**
 * @file TcpDataReceiver.h
 * @created Sep 13, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::daq::TcpDataReceiver
 * Created on: Sep 12, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * (c) Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * (c) Grand Accelerateur National d'Ions Lourds (Ganil)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Frederic Saillant (saillant@ganil.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the MDaq framework.
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

#ifndef mdaq_daq_TcpDataReceiver_h_INCLUDED
#define mdaq_daq_TcpDataReceiver_h_INCLUDED

#include "mdaq/daq/TcpDataTransmit.h"
#include "mdaq/utl/CmdLineArgs.h"
#include "DataReceiver.h"

#include ICE_HEADER(Config.h)
#ifdef ICEE_INT_VERSION
#include ICE_HEADER(Thread.h)
#else
#include <IceUtil/Thread.h>
#endif

#include <vector>

namespace mdaq
{
namespace daq
{

class TcpDataReceiver : virtual public DataReceiver, public TcpDataTransmit
{
	friend class RcveThread;

private:
	struct RcveThread : virtual public ::IceUtil::Thread
	{
		RcveThread(TcpDataReceiver& receiver);
		virtual void run(); ///< Overrides method from ::Ice::Thread
		TcpDataReceiver& receiver;
		std::vector<uint8_t> data;
		volatile bool shouldStop;
		mdaq::utl::CmdLineArgs blips;
		Milliseconds connectionWait_ms;
		Milliseconds receptionWait_ms;
	};
	typedef ::IceUtil::Handle < RcveThread > RcveThreadPtr;

public:
	TcpDataReceiver(const ::utl::net::SocketAddress& sockAddr, size_t rcveBufSize = ~0u);
	virtual ~TcpDataReceiver();

	virtual void start(); // Overrides method from DataReceiver
	virtual void askForStop(); // Overrides method from DataReceiver
	virtual void waitForStop(); // Overrides method from DataReceiver

	Milliseconds connectionWait_ms() const { return rcveThread->connectionWait_ms; }
	Milliseconds receptionWait_ms() const { return rcveThread->receptionWait_ms; }
	void set_connectionWait_ms(Milliseconds _connectionWait_ms) { rcveThread->connectionWait_ms = _connectionWait_ms; }
	void set_receptionWait_ms(Milliseconds _receptionWait_ms) { rcveThread->receptionWait_ms = _receptionWait_ms; }

private:
	int rcveBufSize_;
	const ::utl::net::SocketAddress sockAddr;
	RcveThreadPtr rcveThread;
};

} // namespace daq
} // namespace mdaq

#endif /* mdaq_daq_TcpDataReceiver_h_INCLUDED */
