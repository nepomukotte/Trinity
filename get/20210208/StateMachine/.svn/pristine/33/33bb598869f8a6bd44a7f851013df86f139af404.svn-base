// =============================================================================
// Semaphore.cpp
// -----------------------------------------------------------------------------
// class dhsm::Semaphore
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "Semaphore.h"

namespace dhsm
{
// =============================================================================
// dhsm::Semaphore class methods body
// -----------------------------------------------------------------------------

/**
 * Constructor: initializes the semaphore counter.
 */
Semaphore::Semaphore(int count) : monitor(), count(count)
{
}

/**
 * Destructor.
 */
Semaphore::~Semaphore()
{
}

/**
 * Increments the semaphore counter.
 */
void Semaphore::release()
{
	IceUtil::Monitor<IceUtil::Mutex>::Lock lock(monitor);
	++count;
	monitor.notify();
}

/**
 * Waits for the semaphore counter to be above 0 and then decrements it.
 */
void Semaphore::acquire()
{
	IceUtil::Monitor<IceUtil::Mutex>::Lock lock(monitor);
	while(count <= 0)
	{
		monitor.wait();
	}
}

}
