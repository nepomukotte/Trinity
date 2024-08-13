// =============================================================================
// Engine.cpp
// -----------------------------------------------------------------------------
// class dhsm::Engine
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "Engine.h"
#include "State.h"
#include "Exception.h"
#include <utl/Logging.h>
#include <boost/utility.hpp>
#include <boost/version.hpp>

namespace dhsm
{

// =============================================================================
// dhsm::Engine::StateDataDisposer class methods body
// -----------------------------------------------------------------------------

class Engine::StateDataDisposer : private boost::noncopyable
{
public:
	/**
	 * Clears the data associated with the current state of the specified state machine.
	 */
	StateDataDisposer(StateMachine& sm)
		: sm_(sm)
	{
	}
	~StateDataDisposer()
	{
		try
		{
			Engine::clearStateData(sm_);
		}
		catch(const std::exception& e)
		{
			LOG_ERROR() << "Exception in destructor, state-machine: " << sm_.label() << " message: " << e.what();
		}
	}

private:
	StateMachine& sm_;
};

// =============================================================================
// dhsm::Engine class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
Engine::Engine()
	: stateMachine_()
{
}

/**
 * Destructor. Resets the asssociated state machine, and then destroys it.
 */
Engine::~Engine()
{
	releaseStateMachine();
}

/**
 * @fn StateMachine* Engine::stateMachine()
 * Returns the state machine controlled by this engine.
 */

/**
 * Sets the state machine controlled by this engine, and takes its ownership.
 */
void Engine::setStateMachine(StateMachinePtr sm)
{
	stateMachine_.reset(sm.release());
}

/**
 * Release ownership on the state machine associated with this engine and resets it.
 */
StateMachinePtr Engine::releaseStateMachine()
{
	if(stateMachine_.get() != 0)
	{
		reset(*stateMachine_);
	}
	return stateMachine_;
}

/**
 * Sets the top-level state data.
 * Contrarily to other state data object, which are automatically instanciated,
 * this one can be set by the library-user, giving him complete freedom concerning
 * its initialization. Consequentely, it make it easier for the state machine implementation
 * to communicate with outside objects.
 */
void Engine::setTopContext(StateDataPtr ctx)
{
	topContext_ = ctx;
}

/**
 * @fn State& Engine::currentState(StateMachine& sm)
 * Returns the current state of the specified state machine.
 */

/**
 * Returns the current state at the specified nesting level.
 * @param level the nesting level limit (0 is the outer-most state machine).
 * By default with a level of INT_MAX this method will return the inner-most
 * current state.
 */
State& Engine::currentState(int level)
{
	State* state = stateMachine_->currentState_;
	while(state->isComposite() and level > 0)
	{
		state = state->stateMachine().currentState_;
		--level;
	}
	return *state;
}

/**
 * @fn StateMachine& currentStateMachine(int level)
 * Returns the current state at the specified nesting level.
 * @param level the nesting level limit (0 is the outer-most state machine).
 * By default with a level of INT_MAX this method will return the inner-most
 * current state machine.
 */

/**
 * Returns an active state-machine matching the specified label or 0 if none were found.
 */
StateMachine* Engine::searchCurrentStateMachine(const std::string& smLabel)
{
	StateMachine* sm = stateMachine_.get();
	while(sm != 0 and sm->label() != smLabel)
	{
		sm = currentState(*sm).getStateMachine();
	}
	return sm;
}

StateMachine& Engine::findCurrentStateMachine(const std::string& smLabel)
{
	StateMachine* sm = searchCurrentStateMachine(smLabel);
	if(sm == 0)
	{
		throw StaticError(StaticError::NoSuchStateMachine, smLabel);
	}
	return *sm;
}

/**
 * @fn bool Engine::isInitiated(StateMachine& sm) const
 * Returns true if StateMachine is on, i.e. the state machine's current state is not
 * OFF_STATE.
 */


/**
 * @fn bool Engine::inErrorState(StateMachine& sm) const
 * Returns whether the specified state machine is in its error state.
 */

/**
 * Returns whether the current state machine is in its error state.
 */
bool Engine::inErrorState()
{
	return inErrorState(currentStateMachine());
}

/**
 * Starts the top level state machine.
 * It triggers the crossing of the transition from the initial pseudo state.
 *
 * Nothing happens if the current state machine is already active.
 */
Engine::ProcessingResult Engine::start()
{
	SLOG_DEBUG_BEGIN(msg) << "START END";
	LOG_DEBUG() << "START(" << label_ << ')';
	return doStart();
}

/**
 * Resets and initiates the specified state machine.
 */
Engine::ProcessingResult Engine::reinitiate(StateMachine& sm)
{
	SLOG_DEBUG_BEGIN(msg) << "REINIT END";
	LOG_DEBUG() << "REINIT BEGIN(" << label_ << ',' << sm.path() << ')';
	return doReinitiate(sm);
}

/**
 * Halts the current state machine by putting it in error state.
 * If it throws an std::exception it will put the parent state-machine in error.
 */
Engine::ProcessingResult Engine::halt(StateMachine& sm)
{
	if(not inErrorState(sm))
	{
		SLOG_DEBUG_BEGIN(msg) << "HALT END";
		LOG_DEBUG() << "HALT BEGIN(" << label_ << ',' << sm.path() << ')';
		return doHalt(sm);
	}
	else // we do nothing if the state-machine is already in error
	{
		return Ignored;
	}
}

/**
 * Processes the specified event.
 */
Engine::ProcessingResult Engine::step(const Event& event)
{
	SLOG_DEBUG_BEGIN(endMsg) << "STEP END: " << event.label();
	LOG_DEBUG() << "STEP BEGIN(" << label_ << ',' << event.label() << ',' << stateMachine_->currentState_->path() << ')';
//	checkSupported(event);
	ProcessingResult result;
	try
	{
		result = process(*stateMachine_, event);
	}
	catch(const RecoveredError& e)
	{
		throw TransitionError(e.getMessage());
	}
	catch(const CriticalError& e)
	{
		halt(*stateMachine_); // should not throw, or it would mean entering the top state-machine error state fails.
		throw;
	}
	SLOG_DEBUG(endMsg) << " result: " << result;
	return result;
}


/**
 * Starts the top level state machine.
 * It triggers the crossing of the transition from the initial pseudo state.
 *
 * Nothing happens if the current state machine is already active.
 */
Engine::ProcessingResult Engine::doStart()
{
	initSupportedEvents();
	return initiate(*stateMachine_);
}

/**
 * Starts the specified state machine.
 */
Engine::ProcessingResult Engine::doReinitiate(StateMachine& sm)
{
	reset(sm);
	return initiate(sm);
}

/**
 * Halts the current state machine by putting it in error state.
 * If it throws an std::exception it will put the parent state-machine in error.
 */
Engine::ProcessingResult Engine::doHalt(StateMachine& sm)
{
	onCurrentStateFailsafeExit(sm);
	sm.currentState_ = &sm.errState();
	onCurrentStateEnter(sm, Event());
	return Processed;
}

/**
 * Makes the state machine cross the specified transition.
 */
Engine::ProcessingResult
Engine::crossTransition(StateMachine& sm, const Event& event, Transition& transition)
{
	MLOG_DEBUG("TRANS: " << transition)
	assert(transition.sourceState() == sm.currentState_);

	onCurrentStateExit(sm);

	// execute the transition actions
	onTransitionTrigger(transition, sm.parentCurrentStateData(), event);
	sm.currentState_ = transition.targetState();

	return onCurrentStateEnter(sm, event);
}

/**
 * Starts the specified state machine.
 */
Engine::ProcessingResult Engine::initiate(StateMachine& sm)
{
	if(not isInitiated(sm))
	{
		SLOG_DEBUG_BEGIN(msg) << "INIT END";
		LOG_DEBUG() << "INIT BEGIN(" << label_ << ',' << sm.path() << ')';

		Event event(StateMachine::INITIATE_EVENT_LABEL);
		Transition* transition = sm.currentState_->getTransition(event);
		if(transition != 0)
		{
			return Engine::crossTransition(sm, event, *transition);
		}
		else
		{
			throw StaticError(StaticError::MissingInitialTransition, sm.path());
		}
	}
	else
	{
		return Ignored;
	}
}

/**
 * Recursively processes an event starting from the specified state machine.
 */
Engine::ProcessingResult
Engine::process(StateMachine& sm, const Event& event)
{
	try
	{
		Transition* transition = sm.currentState_->getTransition(event);
		if(transition != 0)
		{
			return crossTransition(sm, event, *transition);
		}
		else
		{
			State* currentState = sm.currentState_;
			if(currentState->isComposite())
			{
				return process(currentState->stateMachine(), event); // recurse one level deeper
			}
			else
			{
				return Ignored;
			}
		}
	}
	catch(const EngineTraversalException& e)
	{
		// these exceptions have to pass through, although they inherit from std::exception
		throw;
	}
	catch(const std::exception& e)
	{
		LOG_DEBUG() << "Exception in state " << sm.currentState_->path() << ": " << e.what();
		halt(sm); // can throw to propagate error to the parent state machine

		// Here we have to convey the exception through the RecoveredError class
		// (which inherits from EngineTraversalException) to allow it to go
		// through the handlers down to the step() method without being caught again.
		// Otherwise the std::exception would be caught at each level, halting
		// every state-machine on the way. So we would lose the ability to stop
		// the faulty sub-state-machine only.
		throw RecoveredError(e.what());
	}
}

/**
 * Resets the specified (sub) state machine.
 * The state-machine is then in OFF state.
 */
Engine::ProcessingResult Engine::reset(StateMachine& sm)
{
	if(isInitiated(sm))
	{
		LOG_DEBUG() << "RESET(" << label_ << ',' << sm.path() << ')';
		onCurrentStateFailsafeExit(sm);
		onStateMachineReset(sm, sm.parentCurrentStateData());
		onCurrentStateEnter(sm, Event());
		return Processed;
	}
	else
	{
		return Ignored;
	}
}

/**
 * This method is called whenever a state machine enters a state.
 * It executes the list of actions associated with entering this State.
 */
Engine::ProcessingResult
Engine::onCurrentStateEnter(StateMachine& sm, const Event& event)
{
	createStateData(sm);
	sm.currentState_->enterActions().execute(event, sm.currentStateData_.get());
	sm.currentStateIsReady_ = true;
	if(sm.currentState_->isComposite())
	{
		return initiate(sm.currentState_->stateMachine());
	}
	return Processed;
}

/**
 * This method is called whenever a state machine exits a state.
 * It executes the list of actions associated with exiting this State.
 */
void Engine::onCurrentStateExit(StateMachine& sm)
{
	StateDataDisposer d(sm); // enable exception-safe disposal of the state-data associated with the current state.

	if(sm.currentState_->isComposite())
	{
		StateMachine& innerSm = sm.currentState_->stateMachine();
		if(innerSm.historyMode() == StateMachine::NO_HISTORY)
		{
			reset(innerSm);
		}
	}

	if(sm.currentStateIsReady_)
	{
		sm.currentStateIsReady_ = false;
		// @TODO: is it correct to create an empty event ???
		sm.currentState_->exitActions().execute(Event(), sm.currentStateData_.get());
	}
}

/**
 * This method is called whenever a state machine exits a state but requires
 * that only critical exceptions can be raised.
 */
void Engine::onCurrentStateFailsafeExit(StateMachine& sm)
{
	try
	{
		onCurrentStateExit(sm);
	}
	catch(const CriticalError& e)
	{
		throw;
	}
	catch(const std::exception& e)
	{
		// Non-critical exceptions raised when exiting a state while resetting or halting a
		// state-machine should not be propagated, because the reset is precisely supposed
		// fix any problem occurring in the state-machine and its states (and sub-statemachines).
		// If the problem is really compromising the upper state-machine, the user code
		// should detect the problem when entering the error state or when resetting the
		// state-machine and throw an exception to put the upper state-machine in error, or
		// throw a CriticalError.
		LOG_WARN() << "Exception leaving state " << sm.currentState_->path() << ": " << e.what();
	}
}

/**
 * This method is called whenever a transition is triggered.
 * It calls successively all the Actions associated with the Transition.
 */
void Engine::onTransitionTrigger(Transition& transition, StateData* closureData, const Event& event)
{
	transition.actions().execute(event, closureData);
}

/**
 * This method is called whenever a state machine is reset.
 * It executes the actions associated with the reset of the state machine.
 */
void Engine::onStateMachineReset(StateMachine& sm, StateData* closureData)
{
	sm.resetActions().execute(Event(), closureData);
	sm.currentState_ = sm.OFF_STATE;
}

/**
 * Clears the data associated with the current state of the specified state machine.
 */
void Engine::clearStateData(StateMachine& sm)
{
	sm.currentStateData_.reset();
}

/**
 * Creates the data associated with the current state of the specified state machine.
 */
void Engine::createStateData(StateMachine& sm)
{
	StateData* context = sm.parentState() != 0 ?
			sm.parentCurrentStateData() :
			topContext_.get();
	sm.currentStateData_.reset(sm.currentState_->createStateData(context, this).release());
}

/**
 * Initializes the list of events supported by the state machine.
 */
void Engine::initSupportedEvents()
{
#ifndef NDEBUG
	supportedEvents.clear();
	listSupportedEvents(*stateMachine_);
#endif
}

/**
 * Checks if the specified event is supported by the state machine.
 */
void Engine::checkSupported(const Event& event)
{
#ifndef NDEBUG
	if(supportedEvents.find(event) == supportedEvents.end())
	{
		throw StaticError(StaticError::UnsupportedEvent, event.label());
	}
#else
	(void) event;
#endif
}

#ifndef NDEBUG
/**
 * Recursive method to build the list of events supported by the state machine
 * (for checking purpose).
 */
void Engine::listSupportedEvents(StateMachine& sm)
{
	for(StateMachine::StateMap::const_iterator currentState = sm.stateMap_.begin(), endState = sm.stateMap_.end();
		currentState != endState; ++currentState)
	{
		State& state = *currentState->second;
		for(State::TransitionMap::const_iterator currentTransition = state.transitionMap().begin(),
			endTransition = state.transitionMap().end();
			currentTransition != endTransition; ++currentTransition)
		{
			const Event& event = currentTransition->first;
			supportedEvents.insert(event);
		}
		if(state.isComposite())
		{
			listSupportedEvents(*state.getStateMachine());
		}
	}
}
#endif

}
