/**
 * @file FdtDataReceiver.h
 * @date March 25th, 2013
 * @author sizun
 *
 * @note
 * SVN tag: $Id$
 * Contributors:
 *   - Shebli Anvar (shebli.anvar@cea.fr)
 *   - Frederic Saillant (saillant@ganil.fr)
 *   - Patrick Sizun (patrick.sizun@cea.fr)
 * @note
 * @note This file is part of the GetBench software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * @copyright © Grand Accelerateur National d'Ions Lourds (Ganil)
 *
 * @par FREE SOFTWARE LICENCING
 *   This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 *   software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 *   CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 *   and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 *   and the software's author, the holder of the economic rights, and the successive licensors have only limited
 *   liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 *   and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 *   mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 *   experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 *   software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 *   to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 *   you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 * @par COMMERCIAL SOFTWARE LICENCING
 *   You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 *   negotiate a specific contract with a legal representative of CEA.
 *
 */

#ifndef mdaq_daq_FdtDataReceiver_h_INCLUDED
#define mdaq_daq_FdtDataReceiver_h_INCLUDED

#include "mdaq/daq/TcpDataTransmit.h"
#include "mdaq/utl/CmdLineArgs.h"
#include "mdaq/daq/DataReceiver.h"
#include <Ice/Ice.h>
#include <IceUtil/Thread.h>
#include <vector>

namespace get
{
namespace daq
{

/**
 * A class to received data sent using the FDT protocol.
 * The receiver filters out heart beats, before transferring data to the data processor.
 */
class FdtDataReceiver : virtual public mdaq::daq::DataReceiver, public mdaq::daq::TcpDataTransmit
{
	friend class RcveThread;

private:
	struct RcveThread : virtual public ::IceUtil::Thread
	{
		RcveThread(FdtDataReceiver& receiver);
		virtual void run(); ///< Overrides method from Ice::Thread
		static size_t decodeFrameSize_B(const uint32_t & chunk);
		FdtDataReceiver& receiver;
		std::vector<uint8_t> data;
		uint32_t word;
		volatile bool shouldStop;
		mdaq::utl::CmdLineArgs blips;
		bool frameStarted;
		size_t remainingDataSize_B;
		Milliseconds connectionWait_ms;
		Milliseconds heartBeatReceptionWait_ms;
		Milliseconds frameReceptionWait_ms;
	};
	typedef ::IceUtil::Handle < RcveThread > RcveThreadPtr;

public:
	FdtDataReceiver(const ::utl::net::SocketAddress& sockAddr, size_t rcveBufSize = ~0u);
	virtual ~FdtDataReceiver();

	virtual void start(); // Overrides method from DataReceiver
	virtual void askForStop(); // Overrides method from DataReceiver
	virtual void waitForStop(); // Overrides method from DataReceiver

private:
	int rcveBufSize_;
	const ::utl::net::SocketAddress sockAddr;
	RcveThreadPtr rcveThread;
};

} // namespace daq
} // namespace get

#endif /* get_daq_FdtDataReceiver_h_INCLUDED */
