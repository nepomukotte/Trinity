// =============================================================================
// Transition.cpp
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

#include "Transition.h"
#include "Event.h"
#include "State.h"
#include "utl/Logging.h"
#include <sstream>

namespace dhsm
{
// =============================================================================
// dhsm::Transition class methods body
// -----------------------------------------------------------------------------

/**
 * To be well defined, any Transition must have a triggering Event, a source
 * State in which it can be triggered and a target State that represents the
 * resulting State of the state machine after a successful Transition.
 * @param triggerEvent The Event that triggers the Transition.
 * @param sourceState The State in which the Transition can be triggered.
 * @param targetState The new State of the state machine if the Transition is
 * successful (does not fail).
 * @note The source and target States can be identical.
 */
Transition::Transition(
	const Event& triggerEvent,
	State& sourceState,
	State& targetState
)
	: triggerEvent_(triggerEvent),
  	sourceState_(&sourceState),
  	targetState_(&targetState),
  	actions_()
{
#ifdef UTL_LOGGING_ENABLED
	std::ostringstream stream;
	stream << "[*]:" << *this;
	actions_.label(stream.str());
#endif
}

#ifdef UTL_LOGGING_ENABLED
void Transition::format(std::ostream& stream) const
{
	stream << '[' << triggerEvent_ << ':' << *sourceState_ << "->" << *targetState_ << ']';
}
#endif

} // namespace dhsm
