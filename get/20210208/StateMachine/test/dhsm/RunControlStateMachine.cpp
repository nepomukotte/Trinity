// =============================================================================
// RunControlStateMachine.cpp
// -----------------------------------------------------------------------------
// class dhsm::RunControlStateMachine
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#define USE_LOG

#include "RunControlStateMachine.h"
#include "dhsm/Exception.h"
#include "utl/Logging.h"
#include <unistd.h>

#if defined(_WIN32)
#include <windows.h>

void sleep(unsigned nsecs)
{
	Sleep(nsecs * 1000);
}
#endif

namespace dhsm
{

// =============================================================================
// dhsm::MasterEngine class methods body
// -----------------------------------------------------------------------------

RunControlStateMachine::RunControlStateMachine(Behaviour behaviour)
	: StateMachine("SmRunControl"), behaviour(behaviour)
{
	state("Off").enterActions().add(*this, &RunControlStateMachine::onTurnedOff);
	state("On").enterActions().add(*this, &RunControlStateMachine::onTurnedOn);
	initial("Off");
	transition("EvTurnOn", "Off", "On");
	transition("EvTurnOff", "On", "Off");

	StateMachine& smOn = state("On").stateMachine("SmOn");
	smOn.state("Idle").enterActions().add(*this, &RunControlStateMachine::onIdling);
	smOn.state("Active").enterActions().add(*this, &RunControlStateMachine::onActivated);
	smOn.initial("Idle");
	smOn.transition("EvStart", "Idle", "Active");
	smOn.transition("EvStop", "Active", "Idle");
	smOn.errState().enterActions().add(*this, &RunControlStateMachine::onError);

	StateMachine& smActive = smOn.state("Active").stateMachine("SmActive");
	smActive.state("Running").enterActions().add(*this, &RunControlStateMachine::onRunning);
	smActive.state("Paused").enterActions().add(*this, &RunControlStateMachine::onPaused);
	smActive.initial("Running");
	smActive.transition("EvPause", "Running", "Paused");
	smActive.transition("EvResume", "Paused", "Running");
}

RunControlStateMachine::~RunControlStateMachine()
{
}

void RunControlStateMachine::onTurnedOn()
{
	LOG_INFO() << "Turned on";
}

void RunControlStateMachine::onTurnedOff()
{
	LOG_INFO() << "Turned off";
}

void RunControlStateMachine::onIdling()
{
	LOG_INFO() << "Idling";
}

void RunControlStateMachine::onActivated()
{
	LOG_INFO() << "Activated";
}

void RunControlStateMachine::onRunning()
{
	if(behaviour == WillFail)
	{
		LOG_INFO() << "Error. Throwing exception...";
		throw std::runtime_error("can't run");
	}
	else if(behaviour == WillCrash)
	{
		LOG_INFO() << "Critical failure. Leaving...";
		exit(EXIT_FAILURE);
	}
	else if(behaviour == WillHang)
	{
		LOG_INFO() << "Hanging for 10 sec...";
		sleep(10);
	}
	else if(behaviour == WillThrowCritical)
	{
		LOG_INFO() << "Throwing a critical error.";
		throw MyCriticalError();
	}
	else
	{
		LOG_INFO() << "Running";
	}
}

void RunControlStateMachine::onPaused()
{
	LOG_INFO() << "Paused";
}

void RunControlStateMachine::onError()
{
	LOG_INFO() << "Faulty";
}

// =============================================================================
// dhsm::RunControlStateMachine::MyCriticalError class methods body
// -----------------------------------------------------------------------------

RunControlStateMachine::MyCriticalError::MyCriticalError()
{
}

RunControlStateMachine::MyCriticalError::~MyCriticalError() throw()
{
}

void RunControlStateMachine::MyCriticalError::writeMessage(std::ostream& stream) const
{
	stream << "just simulating a critical error.";
}

} // namespace dhsm
