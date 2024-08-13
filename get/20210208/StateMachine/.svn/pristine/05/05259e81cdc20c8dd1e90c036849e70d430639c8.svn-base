// =============================================================================
// Transition.h
// -----------------------------------------------------------------------------
// class dhsm::Transition
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_Transition_h_INCLUDED
#define dhsm_Transition_h_INCLUDED

#include "Action.hpp"
#include <ostream>

namespace dhsm
{
// Forward declarations
class Event;

/**
 * This class represents a transition from a source State to a target State,
 * triggered by a specific Event.
 * Any transition object contains a list of Action objects that are executed
 * whenever the transition is triggered.
 * @note The source and target states can be identical.
 */
class Transition
{
public:
	Transition(const Event& triggerEvent, State& sourceState, State& targetState);

	const Event& triggerEvent() const { return triggerEvent_; }
	State* sourceState() const { return sourceState_; }
	State* targetState() const { return targetState_; }
	ActionList& actions() { return actions_; }
	const ActionList& actions() const { return actions_; }

	#ifdef UTL_LOGGING_ENABLED
	void format(std::ostream& stream) const;
	#endif

protected:
	Event triggerEvent_; ///< The Event that triggers the Transition.
	State* sourceState_; ///< The State in which the Transition is activated.
	State* targetState_; ///< The new State of the state machine if the Transition is successful.
	ActionList actions_; ///< The list of Actions that are executed during this Transition.
};

#ifdef UTL_LOGGING_ENABLED
inline std::ostream& operator<<(std::ostream& stream, const Transition& transition)
{
	transition.format(stream);
	return stream;
}
#endif

} // namespace dhsm

#endif // dhsm_Transition_h_INCLUDED
