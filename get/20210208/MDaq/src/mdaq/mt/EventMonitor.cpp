/* =============================================================================
 * EventMonitor.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::mt::EventMonitor
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

#include "EventMonitor.h"
#include <taskLib.h>
#include <stdexcept>

namespace mdaq
{
namespace mt
{
// =============================================================================
// mdaq::mt::EventMonitor class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
EventMonitor::EventMonitor()
	: notifySema(Semaphore::Counting), stoppedSema(Semaphore::Binary, 1)
{
}

EventMonitor::EventMonitor(ProducerType producerType, int delay_us)
	: notifySema(Semaphore::Counting), stoppedSema(Semaphore::Binary, 1)
{
	consumer.reset(new EventConsumer(notifySema, stoppedSema));
	switch (producerType)
	{
	case Polling:
		producer.reset(new PollingProducer(notifySema, delay_us));
		break;
	case Interrupt:
		throw std::runtime_error("Invalid InterruptProducer constructor");
	default:
		throw std::runtime_error("Undefined ProducerType");
	}
}

EventMonitor::EventMonitor(ProducerType producerType, uint32_t baseAddress, int code)
	: notifySema(Semaphore::Counting), stoppedSema(Semaphore::Binary, 1)
{
	consumer.reset(new EventConsumer(notifySema, stoppedSema));
	switch (producerType)
	{
	case Interrupt:
		producer.reset(new InterruptProducer(notifySema, baseAddress, code));
		break;
	case Polling:
		throw std::runtime_error("Invalid PollingProducer constructor");
	default:
		throw std::runtime_error("Undefined ProducerType");
	}
}

/**
 * Destructor.
 */
EventMonitor::~EventMonitor()
{
	stop();
}

/**
 */
void EventMonitor::start(int timeOut_ms)
{
	if (stoppedSema.tryTake())
	{
		consumer->reset();
		producer->start();
		if (timeOut_ms < 0)
			consumer->start();
		else
			consumer->start(timeOut_ms);
	}
}

/**
 */
void EventMonitor::stop()
{
	if (not isStopped())
	{
		try
		{
			producer->stop();
		}
		catch (std::runtime_error& e)
		{
			std::cout << "*** Error by stopping the producer" << std::endl;
			std::cout << "*** " << e.what() << std::endl;
		}
		catch(...)
		{
			std::cout << "*** Error by stopping the producer" << std::endl;
		}

		try
		{
			consumer->stop();
		}
		catch (std::runtime_error& e)
		{
			std::cout << "*** Error by stopping the consumer" << std::endl;
			std::cout << "*** " << e.what() << std::endl;
		}
		catch(...)
		{
			std::cout << "*** Error by stopping the consumer" << std::endl;
		}
		waitForStopped(10000);
	}
}

}
}

#endif
