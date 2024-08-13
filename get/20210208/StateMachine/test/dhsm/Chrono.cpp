// =============================================================================
// Chrono.cpp
// -----------------------------------------------------------------------------
// class dhsm::Chrono
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#define USE_LOG

#include "Chrono.h"
#include "dhsm/State.h"
#include "utl/Logging.h"

namespace dhsm
{

// =============================================================================
// dhsm::Chrono class methods body
// -----------------------------------------------------------------------------

const char Chrono::ST_IDLING[] = "idling";
const char Chrono::ST_ACTIVE[] = "active";
const char Chrono::ST_PAUSED[] = "paused";
const char Chrono::EV_TOGGLE[] = "toggle";
const char Chrono::EV_TRIGGR[] = "triggr";
const char Chrono::EV_SUSPND[] = "suspnd";
const char Chrono::EV_RESUME[] = "resume";
const char Chrono::EV_DSABLE[] = "dsable";

Chrono::Chrono()
{
	StateMachinePtr sm(new StateMachine("SmChrono"));
	// Configure state machine
	sm->transition(EV_TOGGLE, ST_IDLING, ST_ACTIVE).actions().add(*this, &Chrono::onStart);
	sm->transition(EV_TOGGLE, ST_ACTIVE, ST_IDLING).actions().add(*this, &Chrono::onStop);
	sm->transition(EV_TOGGLE, ST_PAUSED, ST_IDLING).actions().add(*this, &Chrono::onStop);
	sm->transition(EV_TRIGGR, ST_ACTIVE, ST_ACTIVE);
	sm->transition(EV_TRIGGR, ST_PAUSED, ST_PAUSED);

	sm->state(ST_ACTIVE).exitActions().add(*this, &Chrono::onIncrCount);
	sm->state(ST_ACTIVE).enterActions().add(*this, &Chrono::dump);
	sm->state(ST_ACTIVE).enterActions().add(*this, &Chrono::onLaunch);
	sm->state(ST_PAUSED).exitActions().add(*this, &Chrono::onIncrMiss);
	sm->state(ST_PAUSED).enterActions().add(*this, &Chrono::dump);
	sm->state(ST_PAUSED).enterActions().add(*this, &Chrono::onLaunch);

	// The following are pure Transitions (no Action attached)
	sm->transition(EV_SUSPND, ST_ACTIVE, ST_PAUSED);
	sm->transition(EV_RESUME, ST_PAUSED, ST_ACTIVE);
	sm->transition(EV_DSABLE, ST_IDLING, StateMachine::OFF_STATE_LABEL);

	// Setting state machine's initial State and enable the Engine.
	sm->initial(ST_IDLING, StateMachine::NO_HISTORY);
	setStateMachine(sm);
}

Chrono::~Chrono()
{
	releaseStateMachine();
}

void Chrono::onStart()
{
	LOG_INFO() << "start";
	count = 0;
	miss = 0;
}

void Chrono::onStop()
{
	LOG_INFO() << "stop";
}

void Chrono::onLaunch()
{
	LOG_INFO() << "launch";
	step(EV_TRIGGR);
}

void Chrono::onIncrCount()
{
	LOG_INFO() << "incrCount";
	count++;
}

void Chrono::onIncrMiss()
{
	LOG_INFO() << "incrMiss";
	miss++;
}

void Chrono::dump()
{
	LOG_INFO() << "LOG: count = " << count << "; miss = " << miss;
}

} // namespace dhsm
