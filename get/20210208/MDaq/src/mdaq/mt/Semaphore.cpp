/* =============================================================================
 * Semaphore.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::mt::Semaphore
 * Created on: 28 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
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

#include "Semaphore.h"
#include <stdexcept>

namespace mdaq
{
namespace mt
{
// =============================================================================
// mdaq::mt::Semaphore class methods body
// -----------------------------------------------------------------------------

/**
 * Constructor.
 */
Semaphore::Semaphore(Type type, int initCount)
	: id(0)
{
	switch (type)
	{
	case Binary:
		id = semBCreate(SEM_Q_FIFO, initCount == 0 ? SEM_EMPTY : SEM_FULL);
		break;
	case Counting:
		id = semCCreate(SEM_Q_FIFO, initCount);
		break;
	default:
		throw std::runtime_error("type of Semaphore is not recognized");
		break;
	}
	if(id == 0)
	{
		throw std::runtime_error("Problem in the creation of the Semaphore");
	}
}

/**
 * Destructor.
 */
Semaphore::~Semaphore()
{
	semDelete(id);
}

/**
 * Unlocks the semaphore.
 *  -# When the semaphore is either binary(BinarSema) or counting(CountSema):
 *  The behavior of this call depends on the state of the semaphore counter
 *  before the call:
 *  - The counter's value is 0: if there are one or more threads that are
 *    blocked by a take() call on this semaphore, one thread is
 *    unblocked. Which thread is unblocked depends on the ordering option
 *    of blocked threads: if isPriorityOrder is set, the thread with
 *    highest priority is unblocked, otherwise it is the thread that has
 *    blocked the earliest.
 *    If there are no blocking threads, the counter is incremented to 1.
 *  - The counter's value is greater than 0 but smaller than its maximum
 *    value(maxCount): the counter is incremented by 1.
 *  - The counter has already reached its maximum value(maxCount):
 *    nothing happens.
 *  -# The semaphore is a mutex(MutexSema):
 *  - If the calling thread has not called take() on this semaphore before,
 *    the give() call is invalid.
 *  - Otherwise, the give() call pops the latest take() call. One thread
 *    may lock a mutex many times but it has to make as many give() calls
 *    to unlock the mutex.
 */
void Semaphore::give()
{
	/*if(id->semType != ((UINT8) SEM_TYPE_COUNTING))
	{
		semGive(id);
	}
	else
	{
		taskLock();
		if(((int) id->state.count) < maxCount)
		{
			semGive(id);
		}
		taskUnlock();
	}*/
	semGive(id);
}

/**
 * Returns the value of the semaphore counter (from 0 to MaxCount) at the time
 * of the call.
 */
int Semaphore::count() const
{
	if(id->semType == SEM_TYPE_COUNTING)
		return (int)(id->state.owner == 0);
	else
		return id->state.count;
}

/**
 * Locks the semaphore.
 * If empty, the call will block until semaphore is unlocked.
 * @return always true.
 */
bool Semaphore::take()
{
	return semTake(id, WAIT_FOREVER) != ERROR;
}

/**
 * Locks the semaphore with a timeout in microseconds.
 * If empty, the call will block until semaphore is unlocked or blocking time
 * exceeds specified number of microseconds.
 * @return true if the semaphore succeeded to take a token, false if the timeout elapsed.
 */
bool Semaphore::take_us(int timeout_usec)
{
	return semTake(id, SysClock::microsecsToTicks(timeout_usec)) != ERROR;
}

/**
 * Releases all threads blocked on this semaphore by a take() call.
 */
void Semaphore::flush()
{
	semFlush(id);
}

}
}

#endif

