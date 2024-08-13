// =============================================================================
// ServantCallback.h
// -----------------------------------------------------------------------------
// class dhsm::ServantCallback
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_ServantCallback_h_INCLUDED
#define dhsm_ServantCallback_h_INCLUDED

#include "Semaphore.h"
#include ICE_HEADER(Exception.h)

namespace dhsm
{
// forward declarations
class ServantDescriptor;

/**
 * This class is the internal implementation of the ServantCallback.
 */
class ServantCallbackImpl
{
public:
	ServantCallbackImpl(ServantDescriptor& servant, Semaphore& completionSemaphore);
	void ice_response(int result);
	void ice_exception(const Ice::Exception& ex);

private:
	ServantDescriptor& servant;
	Semaphore& completionSemaphore;
};

/**
 * This class implements the generated callback interfaces for AMI calls of the
 * EngineServant distributed methods.
 */
template<class AMI>
class ServantCallback : public AMI
{
public:
	ServantCallback(ServantDescriptor& servant, Semaphore& completionSemaphore);
	void ice_response(int result);
	void ice_exception(const Ice::Exception& ex);

private:
	ServantCallbackImpl impl;
};

}

#endif // dhsm_ServantCallback_h_INCLUDED
