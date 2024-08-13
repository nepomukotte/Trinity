/* =============================================================================
 * EventMonitor.h
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

#ifndef mdaq_mt_EventMonitor_INCLUDED
#define mdaq_mt_EventMonitor_INCLUDED

#include "Semaphore.h"
#include "SysClock.h"
#include "EventConsumer.h"
#include "EventProducer.h"
#include "InterruptProducer.h"
#include "PollingProducer.h"
#include <memory>
#include <iostream>

namespace mdaq
{
namespace mt
{
/**
 * This base class provides the basic services for waiting for an event (interrupt or polling timer)
 * and for notifying such an event.
 */
class EventMonitor
{
public:
	EventMonitor();
	EventMonitor(ProducerType producerType, int delay_us);
	EventMonitor(ProducerType producerType, uint32_t baseAddress, int code);
	virtual ~EventMonitor();
	void start(int timeOut_ms);
	void stop(void);
	bool isStopped() { return waitForStopped(0); }
	void waitForStopped() { stoppedSema.take(); stoppedSema.give(); }
	bool waitForStopped(int timeout_ms)
	{
		if (stoppedSema.take_us(timeout_ms * 1000))
		{
			stoppedSema.give();
			return true;
		}
		return false;
	}
protected:
	Semaphore notifySema;
	Semaphore stoppedSema;
public:
	std::auto_ptr<EventProducer> producer;
	std::auto_ptr<EventConsumer> consumer;
};

}
}

#endif // mdaq_mt_EventMonitor_INCLUDED
