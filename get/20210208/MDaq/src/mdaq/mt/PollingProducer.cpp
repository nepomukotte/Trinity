/* =============================================================================
 * PollingProducer.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::mt::PollingProducer
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

#include "PollingProducer.h"
#include <stdexcept>

using std::runtime_error;

namespace mdaq
{
namespace mt
{
// =============================================================================
// mdaq::mt::PollingProducer class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
PollingProducer::PollingProducer(Semaphore& notifySema, int delay_us)
	: EventProducer(notifySema), watchdogId(wdCreate()), delay_ticks(SysClock::microsecsToTicks(delay_us))
{
	if (watchdogId == NULL)
		throw runtime_error("There is not enough memory to create the watchdog");
}

/**
 * Destructor.
 * It is maybe not a good idea to put a throw in a destructor. It may produce some memory leaks
 */
PollingProducer::~PollingProducer()
{
	if (wdDelete(watchdogId) == ERROR)
		std::cout << "The watchdog has NOT been deleted properly" << std::endl;
}

void PollingProducer::start()
{
	if ( wdStart(watchdogId, delay_ticks, reinterpret_cast<FUNCPTR>(&dispatchHandler), reinterpret_cast<int>(this) ) == ERROR )
		throw runtime_error("The watchdog has NOT started properly");
}

void PollingProducer::stop()
{
	if (wdCancel(watchdogId) == ERROR)
		throw runtime_error("The watchdog has NOT been canceled properly");
}

bool PollingProducer::conditionTest()
{
	return true;
}

/**
 * Receives polling timer notifications and dispatches them to the handler
 * It then restarts the watchdog to render the polling periodic
 */
int PollingProducer::dispatchHandler(PollingProducer* pollProd)
{
	if (pollProd->conditionTest() == true)
	{
		pollProd->handler();
		pollProd->notify();
	}
	pollProd->start();
	return 0;
}

} //mt
} //mdaq

#endif
