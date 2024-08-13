// =============================================================================
// StateMachine.cpp
// -----------------------------------------------------------------------------
// class dhsm::StateMachine
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "StateMachine.h"
#include "State.h"
#include "Event.h"
#include <utl/Logging.h>
#include <sstream>

namespace dhsm
{
// =============================================================================
// dhsm::StateMachine class methods body
// -----------------------------------------------------------------------------

const char StateMachine::OFF_STATE_LABEL[] = "__OFF__";
const char StateMachine::ERR_STATE_LABEL[] = "__ERR__";
const char StateMachine::INITIATE_EVENT_LABEL[] = "__INIT__";

/**
 * Constructs a state machine, having only its initial pseudo state (called OFF_STATE).
 * @param label it is possible to give a name to this state machine for debugging purposes.
 * @param parent the parent state if this is a nested state machine.
 * Null if this is a top level state machine.
 */
StateMachine::StateMachine(const char* label, State* parent)
	: label_(label),
	OFF_STATE(new State(this, OFF_STATE_LABEL)),
	ERR_STATE(new State(this, ERR_STATE_LABEL)),
	parentState_(parent),
	stateMap_(),
	resetActions_(),
	currentState_(OFF_STATE),
	currentStateIsReady_(true) // we are ready in OFF state
{
	stateMap_.insert(make_pair(std::string(OFF_STATE_LABEL), OFF_STATE));
	stateMap_.insert(make_pair(std::string(ERR_STATE_LABEL), ERR_STATE));
#ifdef UTL_LOGGING_ENABLED
	std::ostringstream stream;
	stream << "[*]:[onReset(" << *this << ")]";
	resetActions_.label(stream.str());
#endif
}

/**
 * Deallocates all objects created on the heap.
 */
StateMachine::~StateMachine()
{
	// De-allocating state map objects (States)
	for(StateMap::iterator it = stateMap_.begin(); it != stateMap_.end() ; ++it)
	{
		delete it->second;
	}
}

StateMachine& StateMachine::topParentStateMachine()
{
	StateMachine* sm = this;
	while(sm->parentState() != 0)
	{
		sm = sm->parentState()->parentStateMachine();
	}
	return *sm;
}

/**
 * Returns the state data associated with the parent composite state.
 */
StateData* StateMachine::parentCurrentStateData()
{
	StateMachine* parentSm = parentStateMachine();
	return (parentSm != 0) ? parentSm->currentStateData_.get() : 0;
}

/**
 * Returns the number of parent state machines above this one.
 */
int StateMachine::nestingLevel() const
{
	int level = 0;
	const StateMachine* sm = this;
	while(sm->parentState() != 0)
	{
		sm = sm->parentState()->parentStateMachine();
		assert(sm != 0); // each state has a parent state machine, so it shouldn't be null
		++level;
	}
	return level;
}

/**
 * Adds a new State to the StateMachine's set of States. If a State with the same
 * label already exists, nothing happens.
 * @return The just created or already existing State.
 */
State& StateMachine::state(const char* stateLabel)
{
	// Creates new state only if not yet created
	StateMap::iterator it = stateMap_.find(stateLabel);
	if(it != stateMap_.end())
	{
		return *it->second;
	}
	else
	{
		State* newState = new State(this, stateLabel);
		stateMap_.insert(make_pair(newState->label(), newState));
		return *newState;
	}
}

/**
 * Returns a reference to the State matching the specified label.
 * @throw StaticError if the state doesn't exist.
 */
State& StateMachine::findState(const std::string& stateLabel)
{
	StateMap::iterator it = stateMap_.find(stateLabel);
	if(it == stateMap_.end())
	{
		throw StaticError(StaticError::NoSuchState, stateLabel);
	}
	return *it->second;
}

/**
 * Returns a pointer to the State matching the specified label, or null if no one does.
 */
State* StateMachine::searchState(const std::string& stateLabel)
{
	StateMap::iterator it = stateMap_.find(stateLabel);
	if(it == stateMap_.end())
	{
		return 0;
	}
	return it->second;
}

/**
 * @fn const State& StateMachine::offState()
 * Returns a reference to the OFF State of the state machine.
 */

/**
 * This method adds an empty Transition triggered by specified Event between
 * specified states.
 * If some states or event do not exist, they will be created.
 * If the Transition already exists, nothing happens.
 */
Transition& StateMachine::transition(
	const char* evt,
	const char* sourceState,
	const char* targetState
)
{
	return transition(Event(evt), state(sourceState), state(targetState));
}

/**
 * This method adds an empty Transition triggered by specified Event between
 * specified states. If Transition already exists, nothing happens.
 */
Transition& StateMachine::transition(
	const Event& evt,
	State& sourceState,
	State& targetState
)
{
	return sourceState.addTransition(evt, targetState);
}

/**
 * Creates the transition between the initial pseudo state and the entry state.
 * @param entryState The state where the state machine begins, if there is no
 * history or if this is the first time it is run.
 * @param history Tells whether this sub state machine should stay in the state it
 * is (HISTORY) when leaving the parent composite state while running,
 * or whether the state should be reset (NO_HISTORY).
 */
Transition& StateMachine::initial(const char* entryState, HistoryMode history)
{
	return initial(state(entryState), history);
}

/**
 * Creates the transition between the initial pseudo state and the entry state.
 * @param entryState The state where the state machine begins, if there is no
 * history or if this is the first time it is run.
 * @param history Tells whether this sub state machine should stay in the state it
 * is (HISTORY) when leaving the parent composite state while running,
 * or whether the state should be reset (NO_HISTORY).
 */
Transition& StateMachine::initial(State& entryState, HistoryMode history)
{
	historyMode_ = history;
	return transition(Event(INITIATE_EVENT_LABEL), *OFF_STATE, entryState);
}

std::string StateMachine::path() const
{
	std::ostringstream str;
	const State* parent = parentState();
	if(parent != 0)
	{
		parent->path(str);
		str << '/';
	}
	str << *this;
	return str.str();
}

#ifdef UTL_LOGGING_ENABLED
/**
 * Printouts info on state machine.
 */
void StateMachine::dump() const
{
	LOG_INFO() << std::string(80,'*');

	for(StateMap::const_iterator it = stateMap_.begin(); it != stateMap_.end(); ++it)
	{
		State& state = *it->second;
		LOG_INFO() << "STATE: " << state;
		for(State::TransitionMap::const_iterator jt = state.transitionMap().begin();
			jt != state.transitionMap().end(); ++jt)
		{
			const Transition& transition = *jt->second;
			LOG_INFO() << "    " << transition.triggerEvent() << ':'
				 << transition.actions()
				 << '(' << transition.actions().count() << ')';
		}
	}
}
#endif

} // namespace dhsm
