/* =============================================================================
 * Timer.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::mt::Timer
 * template class mdaq::mt::Timer
 * Created on: 28 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the Multitasking sub-package of the Mordicus
 * Real-Time Software Development Framework.
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#include "Timer.h"
#include "SigEvent.hpp"
#ifdef VXWORKS
#include <timers.h>
#else
#include <time.h>
#endif

namespace mdaq
{
namespace mt
{

Timer::Timer(int sigNo_) : SigAction(sigNo_), sigEvent(sigNo_, this)
{
	if (timer_create(CLOCK_REALTIME, &sigEvent, &posixTimer) != 0) throw "Failed to create Posix timer";
}

Timer::~Timer()
{
	timer_delete(posixTimer);
}

void Timer::schedule(long _duration_ms, long _period_ms)
{
	bool immediateCall = false;
	if (_period_ms < 0) _period_ms = 0;
	if (_duration_ms <= 0)
	{
		immediateCall = true; // have to call handler immediately
		_duration_ms = _period_ms;
	}

	struct itimerspec timerparams;

	struct timespec& duration(timerparams.it_value);
	duration.tv_sec  = _duration_ms / 1000;
	duration.tv_nsec = (_duration_ms - duration.tv_sec*1000)*1000000;

	struct timespec& period(timerparams.it_interval);
	period.tv_sec  = _period_ms / 1000;
	period.tv_nsec = (_period_ms - period.tv_sec*1000)*1000000;
	timer_settime(posixTimer, 0, &timerparams, NULL);

	if (immediateCall) handler();
}

void Timer::cancel()
{
	struct itimerspec timerparams;
	timerparams.it_value.tv_sec  = 0;
	timerparams.it_value.tv_nsec = 0;
	timer_settime(posixTimer, 0, &timerparams, NULL);
}

long Timer::period_ms()
{
	struct itimerspec timerparams;
	timer_gettime(posixTimer, &timerparams);
	return timerparams.it_interval.tv_sec*1000 + timerparams.it_interval.tv_nsec/1000000;
}

long Timer::timeLeft_ms()
{
	struct itimerspec timerparams;
	timer_gettime(posixTimer, &timerparams);
	return timerparams.it_value.tv_sec*1000 + timerparams.it_value.tv_nsec/1000000;
}

} // namespace mt
} // namespace mdaq
