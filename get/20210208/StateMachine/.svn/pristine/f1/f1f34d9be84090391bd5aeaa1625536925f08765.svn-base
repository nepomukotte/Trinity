// =============================================================================
// Semaphore.h
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

#ifndef dhsm_Semaphore_h_INCLUDED
#define dhsm_Semaphore_h_INCLUDED

#include ICEUTIL_HEADER(Monitor.h)
#include ICEUTIL_HEADER(Mutex.h)

namespace dhsm
{
/**
 * This is a simple implementation of a counting semaphore based on Ice
 * synchronisation primitives.
 */
class Semaphore
{
public:
	Semaphore(int count);
	~Semaphore();

	void release();
	void acquire();

private:
	IceUtil::Monitor<IceUtil::Mutex> monitor;
	int count;
};

}

#endif // dhsm_Semaphore_h_INCLUDED
