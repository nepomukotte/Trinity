/* =============================================================================
 * EventConsumer.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::mt::EventConsumer
 * Created on: 28 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
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

#ifdef VXWORKS

#include "EventConsumer.h"
#include <stdexcept>
#include <iostream>

namespace mdaq
{
namespace mt
{
// =============================================================================
// mdaq::mt::EventConsumer class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
EventConsumer::EventConsumer(Semaphore& notifySema, Semaphore& stoppedSema)
	: notifySema(notifySema), stoppedSema(stoppedSema), timeOut_tk_(WAIT_FOREVER), stop_(true)
{
}

/**
 * Destructor.
 */
EventConsumer::~EventConsumer()
{
}

void EventConsumer::reset()
{
	stop_ = true;
}

/**
 * @fn void EventConsumer::start()
 * Launches the loop monitoring the event (interrupt/polling timer) with no timeout
 */
void EventConsumer::start(void)
{
	timeOut_tk_ = WAIT_FOREVER;
	start_();
}

/**
 * @fn void EventConsumer::start(int timeOut_ms)
 * Launches the loop monitoring the event (interrupt/polling timer) with a timeout given in milliseconds
 */
void EventConsumer::start(int timeOut_ms)
{
	timeOut_tk_ = SysClock::microsecsToTicks(timeOut_ms*1000);
	start_();
}

/**
 * @fn void EventConsumer::stop()
 * Stops monitoring the event.
 */
void EventConsumer::stop(void)
{
	stop_ = true;
	notifySema.flush();
}

/**
 * @fn void EventConsumer::preWaitAction()
 * Optional function called before the wait(). It is meant to be implemented by the user.
 * It is a normal function, which is NOT inside the Interrupt Service Routine (ISR).
 */
void EventConsumer::preWaitAction()
{
}

/**
 * @fn void EventConsumer::postWaitAction()
 * Function called after the wait(). It is meant to be implemented by the user.
 * It is a normal function, which is NOT inside the Interrupt Service Routine (ISR).
 */
void EventConsumer::postWaitAction()
{
}
//__________________________________________________________________________________________________
/**
 * @fn void EventConsumer::timeOutAction()
 * Function called when a timeout occurs while waiting for an event. It is meant to be implemented by the user.
 */
void EventConsumer::timeOutAction()
{
}
//__________________________________________________________________________________________________
/**
 * Wait for an event
 */
bool EventConsumer::wait_()
{
	return notifySema.take_tk(timeOut_tk_);
}


void EventConsumer::start_()
{
	static char anim[] = { '/', '-', '\\', '|'};
	static unsigned waitCount = 0u;

	stop_ = false;
	while (not stop_)
	{
		preWaitAction();
		if (wait_())
		{
			if (not stop_)
			{
				postWaitAction();
			}
		}
		else
		{
			timeOutAction();
			std::cout << std::string(20, '-') << anim[waitCount] << " Waiting for interrupt " << anim[waitCount] << std::string(20, '-') << '\r' << std::flush;
			waitCount = (waitCount+1)%(sizeof(anim)/sizeof(char));
		}
	}
	stoppedSema.give();
}

}
}

#endif
