// =============================================================================
// State.cpp
// -----------------------------------------------------------------------------
// class dhsm::State
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "StateMachine.h"
#include "Transition.h"
#include "State.h"
#include "Event.h"
#include <boost/version.hpp>
#include <sstream>
#include <memory>

namespace dhsm
{
// =============================================================================
// dhsm::State class methods body
// -----------------------------------------------------------------------------

/**
 * A State constructor must be given a label for human readable referencing.
 */
State::State(StateMachine* parentSm, const char* name)
 	: label_(name), transitionMap_(),
  	enterActions_(), exitActions_(),
  	stateMachine_(), parentStateMachine_(parentSm),
  	createStateData_(0)
{
#ifdef UTL_LOGGING_ENABLED
	std::ostringstream stream;
	stream << "[*]:[onEntr(" << *this << ")]";
	enterActions_.label(stream.str());
	stream.str("");
	stream << "[*]:[onExit(" << *this << ")]";
	exitActions_.label(stream.str());
#endif
}

/**
 * Adds a Transition to the Transition map of this State. If a Transition
 * associated to the specified Event already exists, nothing happens, otherwise
 * a new Transition is created.
 * @param triggerEvent The Event that triggers the Transition.
 * @param targetState The next current State of the state machine whenever the
 * Transition is triggered successfully (does not fail).
 * @note the source State of Transition will naturally be this State.
 * @return Already existing or newly created Transition.
 */
Transition& State::addTransition(
	const Event& triggerEvent,
	State& targetState
)
{
	std::auto_ptr<Transition> transition(new Transition(triggerEvent, *this, targetState));
	Transition& result = *transition;
	TransitionMap::iterator it = transitionMap_.find(triggerEvent);
	if(it != transitionMap_.end())
	{
		transitionMap_.replace(it, transition.release());
	}
	else
	{
		transitionMap_.insert(const_cast<Event&>(triggerEvent), transition.release());
	}
	return result;
}

/**
 * Returns Transition pointer if this State is sensitive to specified Event i.e.
 * if the specified Event causes a Transition trigger. If this State is NOT
 * sensitive to the specified Event, the method returns 0.
 */
Transition* State::getTransition(const Event& triggerEvent)
{
	return isSensitive(triggerEvent) ? transitionMap_.find(triggerEvent)->second : 0;
}

/**
 * @fn bool State::isComposite() const
 * Returns whether this state is composite, ie: if it has a sub state machine.
 */

/**
 * Returns the nested state machine, creating one with the specified name if none is defined yet.
 * @throw StaticError if no name is specified and the state is not composite.
 */
StateMachine& State::stateMachine(const char* name)
{
	if(stateMachine_.get() == 0)
	{
		if(name != 0)
		{
			stateMachine_.reset(new StateMachine(name, this));
		}
		else
		{
			throw StaticError(StaticError::NotCompositeState, path());
		}
	}
	return *stateMachine_;
}

/**
 * Sets the nested state machine.
 * It takes the ownership of the specified state machine, and destroys the
 * previous nested state machine.
 */
StateMachine& State::setStateMachine(std::auto_ptr<StateMachine> subSm)
{
	stateMachine_ = subSm;
	stateMachine_->parentState_ = this;
	return *stateMachine_;
}

/**
 * Tests whether the State is sensitive to specified Event.
 * @return true if the specified event causes a Transition.
 */
bool State::isSensitive(const Event& triggerEvent) const
{
	return transitionMap_.count(triggerEvent) > 0;
}

std::string State::path() const
{
	std::ostringstream str;
	path(str);
	return str.str();
}

void State::path(std::ostringstream& s) const
{
	const State* parent = parentStateMachine()->parentState();
	if(parent != 0)
	{
		parent->path(s);
		s << '/';
	}
	s << *this;
}

} // namespace dhsm
