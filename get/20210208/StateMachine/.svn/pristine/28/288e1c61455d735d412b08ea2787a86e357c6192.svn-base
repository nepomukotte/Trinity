// =============================================================================
// Action.cpp
// -----------------------------------------------------------------------------
// class dhsm::AbstractAction
// class dhsm::Action<T>
// class dhsm::ActionList
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "Action.hpp"
#include "utl/Logging.h"

namespace dhsm
{
// =============================================================================
// dhsm::AbstractAction class methods body
// -----------------------------------------------------------------------------

AbstractAction::~AbstractAction()
{
}

/**
 * @fn void AbstractAction::execute(const Event& event, StateData* closureData) const
 * Executes the action.
 */

// =============================================================================
// Methods of class ActionList
// -----------------------------------------------------------------------------

/**
 * Adds an action to the list.
 * This function takes ownership of the action.
 * @return Returns *this to allow chained calls
 */
ActionList& ActionList::add(AbstractAction* action)
{
	actions.push_back(action);
	return *this;
}

/**
 * Executes all actions of the list in the order in which they have been added
 * to the list.
 * @param event the event that triggered the action.
 * @param closureData the state data associated with the current state, and that
 * must be transmitted to closure actions.
 */
void ActionList::execute(const Event& event, StateData* closureData) const
{
	MLOG_DEBUG("EXEC: " << *this)
	for(ActionVector::const_iterator it = actions.begin(); it != actions.end(); ++it)
	{
		const AbstractAction& action = *it;
		MLOG_DEBUG("EXEC: " << action)
		action.execute(event, closureData);
	}
}

/**
 * @fn size_t ActionList::count() const
 * Returns the total number of Actions contained by this ActionList.
 */

#ifdef UTL_LOGGING_ENABLED
/**
 * Printouts the elements of this ActionList.
 */
void ActionList::dump(int indent) const
{
	for(ActionVector::const_iterator kt = actions.begin(); kt != actions.end(); ++kt)
	{
		LOG_INFO() << std::string(indent, ' ') << *kt;
	}
}
#endif

} // namespace dhsm

