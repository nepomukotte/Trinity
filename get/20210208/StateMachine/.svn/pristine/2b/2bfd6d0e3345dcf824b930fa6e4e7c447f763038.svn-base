// =============================================================================
// Camera.cpp
// -----------------------------------------------------------------------------
// struct dhsm::NotShooting
// struct dhsm::Shooting
// struct dhsm::Idle
// struct dhsm::Configuring
// struct dhsm::Focusing
// struct dhsm::Focused
// struct dhsm::Storing
// class dhsm::Camera
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#define USE_LOG

#include "Camera.h"
#include "dhsm/State.h"
#include "utl/Logging.h"

namespace dhsm
{

struct NotShooting : StateData
{
	NotShooting()
	{
		LOG_INFO() << "entering NotShooting state";
	}
	~NotShooting()
	{
		LOG_INFO() << "leaving NotShooting state";
	}
};

struct Shooting : StateData
{
	Shooting()
	{
		LOG_INFO() << "entering Shooting state";
	}
	~Shooting()
	{
		LOG_INFO() << "leaving Shooting state";
	}
	void displayFocused()
	{
		LOG_INFO() << "Displaying picture...";
	}
};

struct Idle : StateData
{
	Idle()
	{
		LOG_INFO() << "entering Idle state";
	}
	~Idle()
	{
		LOG_INFO() << "leaving Idle state";
	}
};

struct Configuring : StateData
{
	Configuring()
	{
		LOG_INFO() << "entering Configuring state";
	}
	~Configuring()
	{
		LOG_INFO() << "leaving Configuring state";
	}
};

struct Focusing : StateData
{
	Focusing()
	{
		LOG_INFO() << "entering Focusing state";
	}
	~Focusing()
	{
		LOG_INFO() << "leaving Focusing state";
	}
	void focus()
	{
		LOG_INFO() << "Focusing...";
	}
};

struct Focused : StateData
{
	Focused()
	{
		LOG_INFO() << "entering Focused state";
	}
	~Focused()
	{
		LOG_INFO() << "leaving Focused state";
	}
};

struct Storing : StateData
{
	Storing()
	{
		LOG_INFO() << "entering Storing state";
	}
	~Storing()
	{
		LOG_INFO() << "leaving Storing state";
	}
};

// =============================================================================
// dhsm::Camera class methods body
// -----------------------------------------------------------------------------

Camera::Camera()
{
	StateMachinePtr sm(new StateMachine("SmCamera"));
	sm->state<NotShooting>("NotShooting");
	sm->state<Shooting>("Shooting");
	sm->initial("NotShooting", StateMachine::HISTORY);
	sm->transition("EvShutterHalf", "NotShooting", "Shooting");
	sm->transition("EvShutterReleased", "Shooting", "NotShooting");

	StateMachine& smNotShooting = sm->state("NotShooting")
		.stateMachine("SmNotShooting");
	smNotShooting.state<Idle>("Idle");
	smNotShooting.state<Configuring>("Configuring");
	smNotShooting.initial("Idle");
	smNotShooting.transition("EvConfig", "Idle", "Configuring");
	smNotShooting.transition("EvConfig", "Configuring", "Idle");

	StateMachine& smShooting = sm->state("Shooting").stateMachine("SmShooting");
	smShooting.state<Focusing>("Focusing").enterActions().add(&Focusing::focus);
	smShooting.state<Focused>("Focused");
	smShooting.state<Storing>("Storing");
	smShooting.initial("Focusing");
	smShooting.transition("EvInFocus", "Focusing", "Focused")
		.actions().add(&Shooting::displayFocused);
	smShooting.transition("EvShutterFull", "Focused", "Storing");

	engine.setStateMachine(sm);
}

void Camera::start()
{
	LOG_INFO() << ">>>> Initiating camera";
	engine.start();
}

void Camera::pressShutterHalf()
{
	LOG_INFO() << ">>>> Pressed shutter half";
	engine.step("EvShutterHalf");
	engine.step("EvInFocus");
}

void Camera::pressShutterFull()
{
	LOG_INFO() << ">>>> Pressed shutter full";
	engine.step("EvShutterFull");
}

void Camera::releaseShutter()
{
	LOG_INFO() << ">>>> Released shutter";
	engine.step("EvShutterReleased");
}

void Camera::pressConfigure()
{
	LOG_INFO() << ">>>> Pressed configure button";
	engine.step("EvConfig");
}

} // namespace dhsm
