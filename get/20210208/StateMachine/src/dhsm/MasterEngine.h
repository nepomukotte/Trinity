// =============================================================================
// MasterEngine.h
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

#ifndef dhsm_MasterEngine_h_INCLUDED
#define dhsm_MasterEngine_h_INCLUDED

#include "Engine.h"
#include "ServantDescriptor.h"
#include <boost/ptr_container/ptr_vector.hpp>
#include ICE_HEADER(ObjectAdapter.h)
#include ICE_HEADER(Communicator.h)

namespace dhsm
{
// typedefs
typedef boost::ptr_vector<ServantDescriptor>	ServantsList;
typedef ServantsList::iterator					ServantIterator;
typedef ServantsList::const_iterator			ServantConstIterator;

/**
 * This class implements the engine of a distributed state machine, running on
 * an array of hosts and/or processes.
 */
class MasterEngine : public Engine
{
public:
	// inner classes
	class IErrorTolerancePolicy;
	class DefaultErrorPolicy;
	class TooManyFailingServants;
	class TooManyLostServants;

	// typedefs
	typedef std::auto_ptr<IErrorTolerancePolicy> ErrorTolerancePolicyPtr;

	// methods
	MasterEngine();
	MasterEngine(const Ice::CommunicatorPtr& ic, bool destroyOnDtor = false);
	MasterEngine(int& argc, char* argv[], bool destroyOnDtor = true);
	virtual ~MasterEngine();

	//@{
	void init(int& argc, char* argv[], bool destroyOnDtor = true);
	void init(const Ice::CommunicatorPtr& ic, bool destroyOnDtor = false);
	Ice::CommunicatorPtr& communicator() { return ic; }
	//@}
	IErrorTolerancePolicy& errorPolicy();
	void setErrorPolicy(ErrorTolerancePolicyPtr& policy) { errorPolicy_ = policy; }
	//@{
	ServantDescriptor& connectServant(const std::string& proxyStr, bool reinitiate = false);
	ServantDescriptor& connectServant(const char* identity, unsigned short port, const char* host = 0, int timeout = -1, bool reinitiate = false);
	void disconnectAllServants();
	size_t servantCount() const { return servants.size(); }
	ServantDescriptor& servant(size_t index) { return servants[index]; }
	ServantDescriptorPtr detachServant(size_t index);
	void attachServant(ServantDescriptorPtr& servant);
	//@}

protected:
	//@{
	virtual ProcessingResult doStart();
	virtual ProcessingResult doReinitiate(StateMachine& sm);
	virtual ProcessingResult doHalt(StateMachine& sm);
	virtual ProcessingResult crossTransition(StateMachine& sm, const Event& event, Transition& transition);
	//@}

	Ice::CommunicatorPtr ic;	///< the Ice communicator needed to call the servants methods.
	bool destroyOnDtor;		///< tells whether to destroy the ICE communicator in the destructor.
	ServantsList servants;	///< the list of servants to whom we send state machine commands.
	ErrorTolerancePolicyPtr errorPolicy_;
};

/**
 * This is an interface used to define policies about error tolerance in a
 * distributed state machine.
 * It should decide when does the status of the servants becomes too serious to
 * continue the state machine execution normally.
 */
class MasterEngine::IErrorTolerancePolicy
{
public:
	IErrorTolerancePolicy();
	virtual ~IErrorTolerancePolicy();

	virtual void reportErrors(ServantsList& servants) = 0;
};

/**
 * This is a simple error tolerance policy which supports a parameterable number
 * faulty servant.
 * The implementation throw an exception if the number of faulty servants
 * exceeds the threshold, whatever the kind of servant or the kind of error.
 */
class MasterEngine::DefaultErrorPolicy : public MasterEngine::IErrorTolerancePolicy
{
public:
	DefaultErrorPolicy(int maxErrorCount = 0, int maxLostCount = 0);
	virtual ~DefaultErrorPolicy();

	virtual void reportErrors(ServantsList& servants);

protected:
	int maxErrorCount;
	int maxLostCount;
};

/**
 * This exception signals that the number of servants stuck in error state is too high.
 */
class MasterEngine::TooManyFailingServants : public Exception
{
public:
	explicit TooManyFailingServants(const std::string& servantMessages);
	~TooManyFailingServants() throw();

	void writeMessage(std::ostream& stream) const;

private:
	std::string servantMessages;
};

/**
 * This exception signals that the number of servants not responding is too high.
 */
class MasterEngine::TooManyLostServants : public CriticalError
{
public:
	explicit TooManyLostServants(const std::string& servantMessages);
	~TooManyLostServants() throw();

	void writeMessage(std::ostream& stream) const;

private:
	std::string servantMessages;
};

}

#endif // dhsm_MasterEngine_h_INCLUDED
