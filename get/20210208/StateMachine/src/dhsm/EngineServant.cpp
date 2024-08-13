// =============================================================================
// EngineServant.cpp
// -----------------------------------------------------------------------------
// class dhsm::EngineServant
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "EngineServant.h"

#define CATCH_AND_THROW_ICE_EXCEP \
	catch(...) \
	{ \
		handleExceptions(); \
		throw; \
	}

namespace dhsm
{
// =============================================================================
// dhsm::EngineServant class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
EngineServant::EngineServant()
{
}

/**
 * Destructor.
 */
EngineServant::~EngineServant()
{
}

/**
 * Sets the engine that will execute received commands.
 * It is mandatory to call this function before using this class to incarnate an
 * ICE object.
 */
void EngineServant::setEngine(EnginePtr& newEngine)
{
	engine_ = newEngine;
}

/**
 * Callback that starts the state machine engine.
 */
int EngineServant::start(const Ice::Current&)
{
	try
	{
		return engine_->start();
	}
	CATCH_AND_THROW_ICE_EXCEP
}

/**
 * Callback that reinitiates the current state machine at the specified nesting level.
 */
int EngineServant::reinitiate(const std::string& stateMachine, const Ice::Current&)
{
	try
	{
		return engine_->reinitiate(engine_->findCurrentStateMachine(stateMachine));
	}
	CATCH_AND_THROW_ICE_EXCEP
}

/**
 * Callback that puts in error the current state machine at the specified nesting level.
 */
int EngineServant::halt(const std::string& stateMachine, const Ice::Current&)
{
	try
	{
		return engine_->halt(engine_->findCurrentStateMachine(stateMachine));
	}
	CATCH_AND_THROW_ICE_EXCEP
}

/**
 * Callback that processes the received events on the engine.
 */
int EngineServant::stepEvent(const std::string& event, const Ice::Current&)
{
	try
	{
		return engine_->step(Event(event.c_str()));
	}
	CATCH_AND_THROW_ICE_EXCEP
}

/**
 * Callback that processes the received messages having a string parameter on the engine.
 */
int EngineServant::stepStringMsg(const std::string& event, const std::string& param, const Ice::Current&)
{
	try
	{
		Message<std::string> msg(event.c_str(), param);
		return engine_->step(msg);
	}
	CATCH_AND_THROW_ICE_EXCEP
}

/**
 * Callback that processes the received messages having an integer parameter on the engine.
 */
int EngineServant::stepIntMsg(const std::string& event, int param, const Ice::Current&)
{
	try
	{
		Message<int> msg(event.c_str(), param);
		return engine_->step(msg);
	}
	CATCH_AND_THROW_ICE_EXCEP
}

void EngineServant::handleExceptions()
{
	try
	{
		throw;
	}
	catch(const CriticalError& e)
	{
		LOG_DEBUG() << "Propagating critical exception through ICE: " << e;
		throw ServantCriticalError(e.getMessage());
	}
	catch(const std::exception& e)
	{
		LOG_DEBUG() << "Propagating recoverable exception through ICE: " << e.what();
		throw ServantError(e.what());
	}
}

}
