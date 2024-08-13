// =============================================================================
// MasterEngine.cpp
// -----------------------------------------------------------------------------
// class dhsm::MasterEngine
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "MasterEngine.h"
#include "EngineServant.h" // warning: this is a Ice generated header
#include "ServantCallback.h"
#include "ServantDescriptor.h"
#include <utl/Logging.h>
#include ICE_HEADER(Connection.h)
#include ICE_HEADER(Properties.h)
#include ICE_HEADER(Initialize.h)
#include ICEUTIL_HEADER(Config.h)
#include <sstream>

namespace dhsm
{
// =============================================================================
// classes declarations
// -----------------------------------------------------------------------------

/**
 * This private class extends the ServantDescriptor to add the ICE proxy which
 * is necessary to send commands to the servant.
 *
 * @author fchateau
 */
class ServantHandle : public ServantDescriptor
{
public:
	class Identity : public ServantDescriptor::Identity
	{
	public:
		Identity(ServantHandle& servant) : servant(servant) {}

	private:
		void write(std::ostream& stream) const;

		ServantHandle& servant;
	};
	ServantHandle(IEngineServantPrx& servantProxy);
	virtual ~ServantHandle();

	const IEngineServantPrx& proxy() const { return proxy_; }
	IEngineServantPrx& proxy() { return proxy_; }

	const Identity& identity() const { return id_; }

private:
	IEngineServantPrx proxy_;	///< the servant proxy
	Identity id_;
};

// =============================================================================
// dhsm::ServantHandle class methods body
// -----------------------------------------------------------------------------

/**
 * Full constructor.
 */
ServantHandle::ServantHandle(IEngineServantPrx& servantProxy)
	: ServantDescriptor(), proxy_(servantProxy), id_(*this)
{
}

/**
 * Destructor.
 */
ServantHandle::~ServantHandle()
{
}

// =============================================================================
// dhsm::ServantHandle::Identity class methods body
// -----------------------------------------------------------------------------

void ServantHandle::Identity::write(std::ostream& stream) const
{
	const std::string& adapterId = servant.proxy()->ice_getAdapterId();
	if(!adapterId.empty())
	{
		stream << adapterId;
	}
	else
	{
#if ICE_INT_VERSION < 30500
		Ice::ConnectionPtr conn = servant.proxy()->ice_getConnection();
		Ice::IPConnectionInfoPtr info = Ice::IPConnectionInfoPtr::dynamicCast(conn->getInfo());
		if(info)
		{
			stream << info->remoteAddress << ':' << info->remotePort;
		}
#else
		stream << "Unknown";
#endif
	}
}

// =============================================================================
// dhsm::Event class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
MasterEngine::MasterEngine()
	: ic(), destroyOnDtor(false), servants()
{
}

MasterEngine::MasterEngine(const Ice::CommunicatorPtr& ic, bool destroyOnDtor)
	: ic(ic), destroyOnDtor(destroyOnDtor), servants()
{
}

MasterEngine::MasterEngine(int& argc, char* argv[], bool destroyOnDtor)
{
	init(argc, argv, destroyOnDtor);
}

/**
 * Destructor.
 */
MasterEngine::~MasterEngine()
{
	if(destroyOnDtor and ic.get() != 0)
	{
		ic->destroy();
	}
}

/**
 * Initializes the network communication layer.
 * One of the init methods must be called before calling anything else.
 * @param argc the number of the command line arguments given to the program.
 * @param argv the list of command line arguments.
 * It will enable the user to define custom ICE properties on the command line.
 * @param destroyOnDtor tells whether the communicator should be destroyed in
 * the destructor or not.
 */
void MasterEngine::init(int& argc, char* argv[], bool destroyOnDtor)
{
#if ICE_INT_VERSION >= 30300
	ic = Ice::initialize(argc, argv);
#else
	Ice::InitializationData initData;
	initData.properties = Ice::createProperties(argc, argv);
	initData.properties->setProperty("Ice.MonitorConnections", "1"); // enable timeout on AMI calls: check timeouts every 1 second.
	ic = Ice::initialize(initData);
#endif
	this->destroyOnDtor = destroyOnDtor;
}

/**
 * Initializes the network communication layer with an ICE communicator.
 * One of the init methods must be called before calling anything else.
 * @param ic the ICE communicator to use for command distribution.
 * @param destroyOnDtor tells whether the communicator should be destroyed in
 * the destructor or not.
 * @note If you are using a version of Ice below 3.3.0, and if you want to have
 * timeouts on distributed calls, you must set the Ice.MonitorConnections to a
 * value in second which will determine the periodicity at which timeouts will
 * be checked.
 * @note The "destroy" method will be called on the ICE communicator in the
 * engine destructor if it has the last reference on it. Otherwise it is up to
 * the framework user to call this method.
 */
void MasterEngine::init(const Ice::CommunicatorPtr& ic, bool destroyOnDtor)
{
	this->ic = ic;
	this->destroyOnDtor = destroyOnDtor;
}

MasterEngine::IErrorTolerancePolicy& MasterEngine::errorPolicy()
{
	if(errorPolicy_.get() == 0) // if no policy is defined, use the default one
	{
		errorPolicy_.reset(new DefaultErrorPolicy);
	}
	return *errorPolicy_;
}

/**
 * Connects to a servant, and adds it as a member of the distributed state machine.
 */
ServantDescriptor& MasterEngine::connectServant(const std::string& proxyStr, bool reinitiate)
{
	Ice::ObjectPrx base = ic->stringToProxy(proxyStr);
	IEngineServantPrx proxy = IEngineServantPrx::checkedCast(base);
	if(not proxy)
	{
		throw std::runtime_error("invalid EngineServant proxy");
	}
	if(isInitiated(*stateMachine_))
	{
		Ice::Int started = proxy->start();
		if(started == Engine::Ignored and reinitiate)
		{
			proxy->reinitiate(stateMachine_->label());
		}
	}
	ServantHandle* servant = new ServantHandle(proxy);
	servants.push_back(servant);
	return *servant;
}

/**
 * Connects to a servant, and adds it as a member of the distributed state machine.
 */
ServantDescriptor& MasterEngine::connectServant(const char* identity, unsigned short port, const char* host, int timeout, bool reinitiate)
{
	std::ostringstream proxyStr;
	proxyStr << identity << ':' << "default -p " << port;
	if(host != 0)
	{
		proxyStr << " -h " << host;
	}
	if(timeout >= 0)
	{
		proxyStr << " -t " << timeout;
	}
	std::cout << proxyStr.str() << std::endl;
	return connectServant(proxyStr.str(), reinitiate);
}

/**
 * Removes all members of the distributed state machine.
 */
void MasterEngine::disconnectAllServants()
{
	servants.clear();
}

/**
 * @fn size_t MasterEngine::servantCount() const
 * Returns the number of servants controlled by this master engine.
 */

/**
 * @fn ServantDescriptor& MasterEngine::servant(size_t index)
 * Returns the descriptor of the servant at the specified index.
 */

/**
 * Detach a servant from the list of servants controlled by this engine.
 */
ServantDescriptorPtr MasterEngine::detachServant(size_t index)
{
	return ServantDescriptorPtr(servants.release(servants.begin() + index).release());
}

/**
 * Attach a servant previously detached from another engine.
 * The attached servant is appended to the end of the list.
 * @note both engine should use the same ICE communicator, otherwise it may
 * cause problems if the communicator of the source engine is destroyed while
 * the destination engine still uses one of its servants.
 */
void MasterEngine::attachServant(ServantDescriptorPtr& servant)
{
	servants.push_back(servant.release());
}

/**
 * Starts the distributed state machine.
 */
MasterEngine::ProcessingResult MasterEngine::doStart()
{
	Semaphore completionSemaphore(1 - servants.size());
	for(ServantIterator it = servants.begin(), end = servants.end(); it != end; ++it)
	{
		ServantHandle& servant = static_cast<ServantHandle&>(*it);
		if(servant.isLost())
		{
			completionSemaphore.release(); // release the semaphore of ignored servants
			continue;
		}
		AMI_IEngineServant_startPtr callback =
			new ServantCallback<AMI_IEngineServant_start>(servant, completionSemaphore);
		servant.proxy()->start_async(callback);
	}
	completionSemaphore.acquire();
	return Engine::doStart();
}

/*
 * (see inherited doc)
 */
MasterEngine::ProcessingResult MasterEngine::doReinitiate(StateMachine& sm)
{
	Semaphore completionSemaphore(1 - servants.size());
	for(ServantIterator it = servants.begin(), end = servants.end(); it != end; ++it)
	{
		ServantHandle& servant = static_cast<ServantHandle&>(*it);
		if(servant.isLost())
		{
			completionSemaphore.release(); // release the semaphore of ignored servants
			continue;
		}
		AMI_IEngineServant_reinitiatePtr callback =
			new ServantCallback<AMI_IEngineServant_reinitiate>(servant, completionSemaphore);
		servant.proxy()->reinitiate_async(callback, sm.label());
	}
	completionSemaphore.acquire();
	return Engine::doReinitiate(sm);
}

/*
 * (see inherited doc)
 */
MasterEngine::ProcessingResult MasterEngine::doHalt(StateMachine& sm)
{
	Semaphore completionSemaphore(1 - servants.size());
	for(ServantIterator it = servants.begin(), end = servants.end(); it != end; ++it)
	{
		ServantHandle& servant = static_cast<ServantHandle&>(*it);
		if(servant.isLost())
		{
			completionSemaphore.release(); // release the semaphore of ignored servants
			continue;
		}
		AMI_IEngineServant_haltPtr callback =
			new ServantCallback<AMI_IEngineServant_halt>(servant, completionSemaphore);
		servant.proxy()->halt_async(callback, sm.label());
	}
	completionSemaphore.acquire();
	return Engine::doHalt(sm);
}

/*
 * (see inherited doc)
 */
MasterEngine::ProcessingResult
MasterEngine::crossTransition(StateMachine& sm, const Event& event, Transition& transition)
{
	typedef const Message<std::string> StringMessage;
	typedef const Message<int> IntMessage;

	StringMessage* strMsg = dynamic_cast<StringMessage*>(&event);
	IntMessage* intMsg = (strMsg == 0) ? dynamic_cast<IntMessage*>(&event) : 0;
	Semaphore completionSemaphore(1 - servants.size());

	for(ServantIterator it = servants.begin(), end = servants.end(); it != end; ++it)
	{
		ServantHandle& servant = static_cast<ServantHandle&>(*it);
		try
		{
			if(not servant.isOk())
			{
				completionSemaphore.release(); // release the semaphore of ignored servants
				continue;
			}
			if(strMsg != 0)
			{
				AMI_IEngineServant_stepStringMsgPtr callback =
					new ServantCallback<AMI_IEngineServant_stepStringMsg>(servant, completionSemaphore);
				servant.proxy()->stepStringMsg_async(callback, event.label(), strMsg->parameter());
			}
			else if(intMsg != 0)
			{
				AMI_IEngineServant_stepIntMsgPtr callback =
					new ServantCallback<AMI_IEngineServant_stepIntMsg>(servant, completionSemaphore);
				servant.proxy()->stepIntMsg_async(callback, event.label(), intMsg->parameter());
			}
			else
			{
				AMI_IEngineServant_stepEventPtr callback =
					new ServantCallback<AMI_IEngineServant_stepEvent>(servant, completionSemaphore);
				servant.proxy()->stepEvent_async(callback, event.label());
			}
		}
		catch(const std::exception& e)
		{
			LOG_ERROR() << "Exception when calling step() on servant proxy. Details\n" << e.what();
			servant.setStatus(ServantDescriptor::Critical);
		}
	}
	completionSemaphore.acquire();

	// the policy is called here because it should throw an exception in case of
	// problem, and this location takes advantage of the master recursive handlers
	// to enable the standard error state mechanism apply to these exceptions as well.
	errorPolicy().reportErrors(servants);
	return Engine::crossTransition(sm, event, transition);
}

// =============================================================================
// dhsm::MasterEngine::IErrorTolerancePolicy class methods body
// -----------------------------------------------------------------------------

MasterEngine::IErrorTolerancePolicy::IErrorTolerancePolicy()
{
}

MasterEngine::IErrorTolerancePolicy::~IErrorTolerancePolicy()
{
}

/**
 * @fn void MasterEngine::IErrorTolerancePolicy::reportErrors(ServantsList& servants)
 * Analyses the error status of each servant and decides to report the error if
 * it is too serious.
 * The framework client must override this method to define custom error policies.
 * This method should use exceptions to report any error considered serious
 * enough for the policy.
 */

// =============================================================================
// dhsm::MasterEngine::DefaultErrorPolicy class methods body
// -----------------------------------------------------------------------------

MasterEngine::DefaultErrorPolicy::DefaultErrorPolicy(int maxErrorCount, int maxLostCount)
	: maxErrorCount(maxErrorCount), maxLostCount(maxLostCount)
{
}

MasterEngine::DefaultErrorPolicy::~DefaultErrorPolicy()
{
}

/*
 * (see inherited doc)
 */
void MasterEngine::DefaultErrorPolicy::reportErrors(ServantsList& servants)
{
	static const char* const statusStr[] = {
		"working", "failing", "critical", "not responding", "unknown state"
	};
	std::ostringstream aggregatedMessages;
	int errorCount = 0, lostCount = 0;
	for(ServantIterator it = servants.begin(), end = servants.end(); it != end; ++it)
	{
		ServantHandle& servant = dynamic_cast<ServantHandle&>(*it);
		if(servant.isBad())
		{
			++errorCount;
		}
		else if(servant.isLost())
		{
			++lostCount;
		}
		if(!servant.isOk())
		{
			aggregatedMessages << std::endl
			 	<< servant.identity() << ' ' << statusStr[servant.status()]
			    << ": " << servant.errorMessage();
		}
	}
	if(lostCount > maxLostCount)
	{
		throw TooManyLostServants(aggregatedMessages.str());
	}
	if(errorCount > maxErrorCount)
	{
		throw TooManyFailingServants(aggregatedMessages.str());
	}
}

// =============================================================================
// dhsm::MasterEngine::TooManyFailingServants class methods body
// -----------------------------------------------------------------------------

MasterEngine::TooManyFailingServants::TooManyFailingServants(const std::string& servantMessages)
	: servantMessages(servantMessages)
{
}

MasterEngine::TooManyFailingServants::~TooManyFailingServants() throw()
{
}

void MasterEngine::TooManyFailingServants::writeMessage(std::ostream& stream) const
{
	stream << "Too many servants are in error state: " << servantMessages;
}

// =============================================================================
// dhsm::MasterEngine::TooManyLostServants class methods body
// -----------------------------------------------------------------------------

MasterEngine::TooManyLostServants::TooManyLostServants(const std::string& servantMessages)
	: servantMessages(servantMessages)
{
}

MasterEngine::TooManyLostServants::~TooManyLostServants() throw()
{
}

void MasterEngine::TooManyLostServants::writeMessage(std::ostream& stream) const
{
	stream << "Too many servants have been lost: " << servantMessages;
}

}
