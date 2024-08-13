/* =============================================================================
 * SysClock.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::mt::SysClock
 * Created on: 2 October 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
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

#include "SysClock.h"
#include <ctime>

#if defined(__APPLE__)
#include <mach/mach_time.h>
#include <cmath>
#endif

namespace mdaq
{
namespace mt
{
/**
 * @fn unsigned long SysClock::getTicks()
 * Returns the number of ticks since the system start.
 */

/**
 * @fn UsecCount SysClock::getResolution_Hz()
 * Returns the resolution of the system clock (in Hertz).
 * This represents the maximum frequency that a timer can achieve on this
 * platform and the granularity of possible frequencies (all effective frequencies
 * fit in the formula: resolution / n).
 */

/**
 * Returns the resolution of the clock (in microseconds).
 */
UsecCount SysClock::getResolution_us()
{
	return utl::divide<utl::RoundDown>(MICROSEC_PER_SEC, getResolution_Hz());
}

template<utl::RoundMode round>
TickCount SysClock::microsecsToTicks(UsecCount microSecs)
{
	return utl::divide<round>(microSecs * getResolution_Hz(), MICROSEC_PER_SEC);
}

template<utl::RoundMode round>
UsecCount SysClock::ticksToMicrosecs(TickCount ticks)
{
	return utl::divide<round>(ticks * MICROSEC_PER_SEC, getResolution_Hz());
}

UsecCount SysClock::alignDuration(UsecCount microSecs)
{
	return ticksToMicrosecs(microsecsToTicks(microSecs));
}

template TickCount SysClock::microsecsToTicks<utl::RoundDown>(UsecCount microSecs);
template TickCount SysClock::microsecsToTicks<utl::RoundNear>(UsecCount microSecs);
template TickCount SysClock::microsecsToTicks<utl::RoundUp  >(UsecCount microSecs);
template UsecCount SysClock::ticksToMicrosecs<utl::RoundDown>(TickCount ticks);
template UsecCount SysClock::ticksToMicrosecs<utl::RoundNear>(TickCount ticks);
template UsecCount SysClock::ticksToMicrosecs<utl::RoundUp  >(TickCount ticks);

#if ((defined(_BSD_SOURCE) || _XOPEN_SOURCE >= 500) && defined(_POSIX_MONOTONIC_CLOCK)) || defined(__APPLE__)
#include <unistd.h>

void SysClock::sleep_us(UsecCount microSecs)
{
	usleep(microSecs);
}
void SysClock::sleep_tk(TickCount ticks)
{
	usleep(ticks);
}

#if defined(__APPLE__)
TickCount SysClock::getTicks()
{
	uint64_t now = mach_absolute_time();
	static double conversion = 0.0;

	if (conversion == 0.0)
	{
		mach_timebase_info_data_t info;
		kern_return_t err = mach_timebase_info(&info);

		// Convert the timebase into microseconds
		if (err == 0)
			conversion = 1e-3 * (double) info.numer / (double) info.denom;
	}

	return static_cast<TickCount>(std::floor(conversion * (double) now));
}
#else
TickCount SysClock::getTicks()
{
	timespec currentTime;
	clock_gettime(CLOCK_MONOTONIC, &currentTime);
	return static_cast<TickCount>(currentTime.tv_sec) * MICROSEC_PER_SEC + (currentTime.tv_nsec / 1000);;
}
#endif

TickCount SysClock::getResolution_Hz()
{
	return MICROSEC_PER_SEC;
}

#elif defined(VXWORKS)

#include <vxWorks.h>
#include <tickLib.h>
#include <sysLib.h>
#include <taskLib.h>

void SysClock::sleep_us(UsecCount microSecs)
{
	taskDelay(microsecsToTicks(microSecs));
}
void SysClock::sleep_tk(TickCount ticks)
{
	taskDelay(ticks);
}
TickCount SysClock::getTicks()
{
	return tickGet();
}
TickCount SysClock::getResolution_Hz()
{
	return sysClkRateGet();
}

#elif defined(_WIN32)

#include <windows.h>

void SysClock::sleep_us(UsecCount microSecs)
{
	Sleep(microSecs / 1000);
}
void SysClock::sleep_tk(TickCount ticks)
{
	Sleep(ticks);
}
TickCount SysClock::getTicks()
{
	return GetTickCount();
}
TickCount SysClock::getResolution_Hz()
{
	return 1000;
}

#endif

}
}

