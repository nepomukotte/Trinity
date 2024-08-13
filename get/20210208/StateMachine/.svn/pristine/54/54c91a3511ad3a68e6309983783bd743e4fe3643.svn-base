// =============================================================================
// ServantCallback.cpp
// -----------------------------------------------------------------------------
// class dhsm::ServantCallback
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "ServantCallback.h"
#include "EngineServant.h"
#include "Engine.h"
#include "ServantDescriptor.h"
#include <utl/Logging.h>
#include ICE_HEADER(LocalException.h)
#include <boost/lexical_cast.hpp>

namespace dhsm
{
// =============================================================================
// dhsm::ServantCallbackImpl class methods body
// -----------------------------------------------------------------------------

/**
 * Constructor.
 */
ServantCallbackImpl::ServantCallbackImpl(ServantDescriptor& servant,
	Semaphore& completionSemaphore)
	: servant(servant), completionSemaphore(completionSemaphore)
{
}

void ServantCallbackImpl::ice_response(int /*result*/)
{
	servant.setStatus(ServantDescriptor::Working);
	completionSemaphore.release();
}

void ServantCallbackImpl::ice_exception(const Ice::Exception& ex)
{
	try
	{
		ex.ice_throw();
    }
	catch(const ServantCriticalError& e)
	{
		LOG_WARN() << "Servant critical error " << e.reason;
		servant.setStatus(ServantDescriptor::Critical);
		servant.setErrorMessage(e.reason);
	}
	catch(const ServantError& e)
	{
		LOG_WARN() << "Servant error " << e.reason;
		servant.setStatus(ServantDescriptor::Failing);
		servant.setErrorMessage(e.reason);
	}
	catch(const Ice::TimeoutException& e)
	{
		LOG_WARN() << "Servant timeout error " << e;
		servant.setStatus(ServantDescriptor::NotResponding);
		servant.setErrorMessage(boost::lexical_cast<std::string>(e));
    }
	catch(const Ice::SocketException& e)
	{
		LOG_WARN() << "Servant socket exception " << e;
		servant.setStatus(ServantDescriptor::NotResponding);
		servant.setErrorMessage(boost::lexical_cast<std::string>(e));
	}
	catch(const Ice::NoEndpointException& e)
	{
		LOG_WARN() << "Servant no endpoint exception " << e;
		servant.setStatus(ServantDescriptor::NotResponding);
		servant.setErrorMessage(boost::lexical_cast<std::string>(e));
	}
	catch(const Ice::Exception& e)
	{
		LOG_WARN() << "Servant exception " << e;
		servant.setStatus(ServantDescriptor::Critical);
		servant.setErrorMessage(boost::lexical_cast<std::string>(e));
	}
	completionSemaphore.release();
}

// =============================================================================
// dhsm::ServantCallback class methods body
// -----------------------------------------------------------------------------

/**
 * Constructor.
 */
template<class AMI>
ServantCallback<AMI>::ServantCallback(ServantDescriptor& servant,
	Semaphore& completionSemaphore)
	: impl(servant, completionSemaphore)
{
}

template<class AMI>
void ServantCallback<AMI>::ice_response(int result)
{
	impl.ice_response(result);
}

template<class AMI>
void ServantCallback<AMI>::ice_exception(const Ice::Exception& ex)
{
	impl.ice_exception(ex);
}

// explicit template instanciations
template class ServantCallback<AMI_IEngineServant_start>;
template class ServantCallback<AMI_IEngineServant_reinitiate>;
template class ServantCallback<AMI_IEngineServant_halt>;
template class ServantCallback<AMI_IEngineServant_stepEvent>;
template class ServantCallback<AMI_IEngineServant_stepStringMsg>;
template class ServantCallback<AMI_IEngineServant_stepIntMsg>;

}
