/**
 * @file NetPerfProcessor.h
 * @created Sep 14, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::daq::NetPerfProcessor
 * Created on: Sep 14, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
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
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef mdaq_daq_NetPerfProcessor_h_INCLUDED
#define mdaq_daq_NetPerfProcessor_h_INCLUDED

#include "mdaq/daq/DataProcessorCore.h"
#include "mdaq/mt/Timer.h"

#include ICE_HEADER(ICE_EDITION.h)
#include ICE_HEADER(Config.h)
#ifdef ICEE_INT_VERSION
#include ICE_HEADER(Time.h)
#else
#include <IceUtil/Time.h>
#endif

namespace mdaq
{
namespace daq
{

class NetPerfProcessor;

class PerfRcveTimer : public mdaq::mt::Timer
{
public:
	PerfRcveTimer(NetPerfProcessor& netPerf) : netPerf(netPerf) {}
	virtual void handler();

private:
	NetPerfProcessor& netPerf;
};

class NetPerfProcessor : public ::mdaq::daq::DataProcessorCore
{
public:
	NetPerfProcessor(int printPeriod_s = 5);
	virtual ~NetPerfProcessor();

	///@{ @override from mdaq::daq::DataProcessorCore
	virtual void resetData();
	virtual void setActive(bool _isActive = true);
	///@}

	virtual void processData(const uint8_t* buffer, uint32_t nBytes); ///< @override from mdaq::daq::DataProcessorCore
	void startTimer();
	void stopTimer();

protected:
	int printPeriod_ms;
	PerfRcveTimer perfRcveTimer;
	double bytesReceived;
	IceUtil::Time startTime;
	static const char* units[];
};

} // namespace daq
} // namespace mdaq

#endif /* mdaq_daq_NetPerfProcessor_h_INCLUDED */
