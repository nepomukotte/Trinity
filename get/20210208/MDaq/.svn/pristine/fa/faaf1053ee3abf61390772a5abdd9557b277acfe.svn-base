/* =============================================================================
 * EventConsumerAdapter.hpp
 * -----------------------------------------------------------------------------
 * template class mdaq::mt::EventConsumerAdapter
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

#ifndef mdaq_mt_EventConsumerAdapter_INCLUDED
#define mdaq_mt_EventConsumerAdapter_INCLUDED

#include "EventConsumer.h"
#include "DefaultForTemplate.h"
#include <stdexcept>

namespace mdaq
{
namespace mt
{

extern DefaultHandlerClass defaultHandlerClass;

template <typename ProcClass = DefaultHandlerClass>
class EventConsumerAdapter : public EventConsumer
{
public:
	typedef void (ProcClass::*Action)();
	typedef void (EventConsumerAdapter::*InternalFunction)();
	EventConsumerAdapter(Semaphore& notifySema, Semaphore& stoppedSema);
	EventConsumerAdapter(Semaphore& notifySema, Semaphore& stoppedSema, ProcClass* actionObject, Action postWaitActionProc);
	EventConsumerAdapter(Semaphore& notifySema, Semaphore& stoppedSema, ProcClass* actionObject, Action preWaitActionProc, Action postWaitActionProc);
	EventConsumerAdapter(Semaphore& notifySema, Semaphore& stoppedSema, ProcClass* actionObject, Action preWaitActionProc, Action postWaitActionProc, Action timeOutActionProc);
	virtual ~EventConsumerAdapter();
	virtual void preWaitAction();
	virtual void postWaitAction();
	virtual void timeOutAction();
private:
	ProcClass* actionObject;
	void runTimeCheck();
	void preWaitActionOn();
	void preWaitActionOff();
	Action preWaitActionProc;
	Action postWaitActionProc;
	Action timeOutActionProc_;
	InternalFunction preWaitActionWrapper;
};

template<typename ProcClass>
EventConsumerAdapter<ProcClass>::EventConsumerAdapter(Semaphore& notifySema, Semaphore& stoppedSema)
	: EventConsumer(notifySema, stoppedSema), actionObject(&defaultHandlerClass),
	  preWaitActionProc(&DefaultHandlerClass::preWaitAction), postWaitActionProc(&DefaultHandlerClass::postWaitAction),
	  timeOutActionProc_(0),
	  preWaitActionWrapper(&EventConsumerAdapter::preWaitActionOn)
{
	runTimeCheck();
}

template<typename ProcClass>
EventConsumerAdapter<ProcClass>::EventConsumerAdapter(Semaphore& notifySema, Semaphore& stoppedSema,
		ProcClass* actionObject, Action postWaitActionProc)
	: EventConsumer(notifySema, stoppedSema), actionObject(actionObject), preWaitActionProc(0), postWaitActionProc(postWaitActionProc),
	  timeOutActionProc_(0),
	  preWaitActionWrapper(&EventConsumerAdapter::preWaitActionOff)
{
	runTimeCheck();
}

template<typename ProcClass>
EventConsumerAdapter<ProcClass>::EventConsumerAdapter(Semaphore& notifySema, Semaphore& stoppedSema,
		ProcClass* actionObject, Action preWaitActionProc, Action postWaitActionProc)
	: EventConsumer(notifySema, stoppedSema), actionObject(actionObject),
	  preWaitActionProc(preWaitActionProc), postWaitActionProc(postWaitActionProc),
	  timeOutActionProc_(0)
{
	if (preWaitActionProc == 0)
		preWaitActionWrapper = &EventConsumerAdapter::preWaitActionOff;
	else
		preWaitActionWrapper = &EventConsumerAdapter::preWaitActionOn;
	runTimeCheck();
}

template<typename ProcClass>
EventConsumerAdapter<ProcClass>::EventConsumerAdapter(Semaphore& notifySema, Semaphore& stoppedSema,
		ProcClass* actionObject, Action preWaitActionProc, Action postWaitActionProc, Action timeOutActionProc)
	: EventConsumer(notifySema, stoppedSema), actionObject(actionObject),
	  preWaitActionProc(preWaitActionProc), postWaitActionProc(postWaitActionProc),
	  timeOutActionProc_(timeOutActionProc)
{
	if (preWaitActionProc == 0)
		preWaitActionWrapper = &EventConsumerAdapter::preWaitActionOff;
	else
		preWaitActionWrapper = &EventConsumerAdapter::preWaitActionOn;
	runTimeCheck();
}


template<typename ProcClass>
EventConsumerAdapter<ProcClass>::~EventConsumerAdapter()
{
}

template<typename ProcClass>
void EventConsumerAdapter<ProcClass>::preWaitAction(void)
{
	(*this.*preWaitActionWrapper)();
}

template<typename ProcClass>
void EventConsumerAdapter<ProcClass>::preWaitActionOn(void)
{
	(actionObject->*preWaitActionProc)();
}

template<typename ProcClass>
void EventConsumerAdapter<ProcClass>::preWaitActionOff(void)
{
}

template<typename ProcClass>
void EventConsumerAdapter<ProcClass>::postWaitAction(void)
{
	(actionObject->*postWaitActionProc)();
}

template<typename ProcClass>
void EventConsumerAdapter<ProcClass>::timeOutAction()
{
	if (timeOutActionProc_)
	{
		(actionObject->*timeOutActionProc_)();
	}
	else
	{
		EventConsumer::timeOutAction();
	}
}

template<typename ProcClass>
void EventConsumerAdapter<ProcClass>::runTimeCheck(void)
{
	if (actionObject == 0)
		throw std::runtime_error("EventConsumerAdapter: actionObject cannot be NULL");
	if (postWaitActionProc == 0)
		throw std::runtime_error("EventConsumerAdapter: postWaitAction cannot be NULL");
}

}
}

#endif // mdaq_mt_EventConsumerAdapter_INCLUDED
