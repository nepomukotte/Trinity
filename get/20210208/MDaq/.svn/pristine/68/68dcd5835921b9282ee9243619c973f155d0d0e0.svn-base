/* =============================================================================
 * PollingProducerAdapter.hpp
 * -----------------------------------------------------------------------------
 * template class mdaq::mt::PollingProducerAdapter
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

#ifndef mdaq_mt_PollingProducerAdapter_INCLUDED
#define mdaq_mt_PollingProducerAdapter_INCLUDED

#include "EventProducer.h"
#include "PollingProducer.h"
#include "DefaultForTemplate.h"
#include <stdexcept>

namespace mdaq
{
namespace mt
{

extern DefaultHandlerClass defaultHandlerClass;

template<typename ProcClass = DefaultHandlerClass>
class PollingProducerAdapter : public PollingProducer
{
public:
	typedef bool (ProcClass::*ConditionTest)();
	typedef void (ProcClass::*Handler)();
	typedef void (PollingProducerAdapter::*InternalFunction)();
	PollingProducerAdapter(Semaphore& notifySema, int delay_us = defaultPollingDelay_us);
	PollingProducerAdapter(Semaphore& notifySema, int delay_us, ProcClass* handlerObject,
			ConditionTest condTestProc);
	PollingProducerAdapter(Semaphore& notifySema, int delay_us, ProcClass* handlerObject,
			ConditionTest condTestProc, Handler handlerProc);
	virtual ~PollingProducerAdapter();
	virtual bool conditionTest();
	virtual void handler();
private:
	void runTimeCheck();
	void handlerOn();
	void handlerOff();
	ProcClass* handlerObject;
	ConditionTest  condTestProc;
	Handler handlerProc;
	InternalFunction handlerWrapper;
};

template<typename ProcClass>
PollingProducerAdapter<ProcClass>::PollingProducerAdapter(Semaphore& notifySema, int delay_us)
	: PollingProducer(notifySema, delay_us), handlerObject(&defaultHandlerClass),
	  condTestProc(&DefaultHandlerClass::isItTrue), handlerProc(&DefaultHandlerClass::handler),
	  handlerWrapper(&PollingProducerAdapter::handlerOn)
{
	runTimeCheck();
}

template<typename ProcClass>
PollingProducerAdapter<ProcClass>::PollingProducerAdapter(Semaphore& notifySema, int delay_us,
		ProcClass* handlerObject, ConditionTest condTestProc)
	: PollingProducer(notifySema, delay_us), handlerObject(handlerObject),
	  condTestProc(condTestProc), handlerProc(0),
	  handlerWrapper(&PollingProducerAdapter::handlerOff)
{
	runTimeCheck();
}

template<typename ProcClass>
PollingProducerAdapter<ProcClass>::PollingProducerAdapter(Semaphore& notifySema, int delay_us,
		ProcClass* handlerObject, ConditionTest condTestProc, Handler handlerProc)
	: PollingProducer(notifySema, delay_us), handlerObject(handlerObject),
	  condTestProc(condTestProc), handlerProc(handlerProc)
{
	if (handlerProc == 0)
		handlerWrapper = &PollingProducerAdapter::handlerOff;
	else
		handlerWrapper = &PollingProducerAdapter::handlerOn;
	runTimeCheck();
}

template<typename ProcClass>
PollingProducerAdapter<ProcClass>::~PollingProducerAdapter()
{
}

template<typename ProcClass>
bool PollingProducerAdapter<ProcClass>::conditionTest(void)
{
	return (handlerObject->*condTestProc)();
}

template<typename ProcClass>
void PollingProducerAdapter<ProcClass>::handler(void)
{
	(*this.*handlerWrapper)();
}

template<typename ProcClass>
void PollingProducerAdapter<ProcClass>::handlerOn(void)
{
	(handlerObject->*handlerProc)();
}

template<typename ProcClass>
void PollingProducerAdapter<ProcClass>::handlerOff(void)
{
}

template<typename ProcClass>
void PollingProducerAdapter<ProcClass>::runTimeCheck(void)
{
	if (handlerObject == 0)
		throw std::runtime_error("PollingProducerAdapter: handlerObject can't be NULL");
	if (condTestProc == 0)
		throw std::runtime_error("PollingProducerAdapter: condTest can't be NULL");
}

}
}

#endif // mdaq_mt_PollingProducerAdapter_INCLUDED
