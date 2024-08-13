/**
 * @file NetPerfProcessor.cpp
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

#include "NetPerfProcessor.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <climits>

using namespace std;

namespace mdaq
{
namespace daq
{

const char* NetPerfProcessor::units[] = {" b/s", " Kb/s", " Mb/s", " Gb/s"};

void PerfRcveTimer::handler()
{
	netPerf.resetData();
}

NetPerfProcessor::NetPerfProcessor(int printPeriod_s)
	: DataProcessorCore(false),
	  printPeriod_ms(1000*printPeriod_s),
	  perfRcveTimer(*this)
{
	cout << '\n' << string(80, '=') << endl;
	cout.precision(1);
	cout << fixed;
	resetData();
}

NetPerfProcessor::~NetPerfProcessor()
{
	cout << '\n' << string(80, '-') << endl;
}

void NetPerfProcessor::startTimer()
{
	perfRcveTimer.schedule(printPeriod_ms, printPeriod_ms);
	resetData();
}

void NetPerfProcessor::stopTimer()
{
	perfRcveTimer.cancel();
	resetData();
}

void NetPerfProcessor::resetData()
{
	if (isActive())
	{
		double perf = bytesReceived / (IceUtil::Time::now() - startTime).toSecondsDouble();
		perf *= CHAR_BIT; // nb of bits per second
		int iunit = 3;
		if      (perf >= 1e15) { perf *= 1e-9;          }
		else if (perf >= 1e+6) { perf *= 1e-6; --iunit; }
		else if (perf >= 1e+3) { perf *= 1e-3; --iunit; }
		else                   {               --iunit; }
		if (perf > 0.0) cout << ">>> " << perf << units[iunit] << " <<<" << endl;
	}

	bytesReceived = 0.;
	startTime = IceUtil::Time::now();
}

void NetPerfProcessor::setActive(bool _isActive)
{
	if (_isActive)
	{
		startTimer();
		DataProcessorCore::setActive(true);
	}
	else
	{
		DataProcessorCore::setActive(false);
		stopTimer();
	}
}

void NetPerfProcessor::processData(const uint8_t* /*buffer*/, uint32_t nBytes)
{
	bytesReceived += nBytes;
}


} // namespace daq
} // namespace mdaq
