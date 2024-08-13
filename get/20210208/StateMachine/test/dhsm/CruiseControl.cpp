// =============================================================================
// CruiseControl.h
// -----------------------------------------------------------------------------
// class dhsm::HardwareFailure
// struct dhsm::CruiseControlData
// struct dhsm::CruiseData
// class dhsm::CruiseControl
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#define USE_LOG

#include "CruiseControl.h"
#include "dhsm/State.h"
#include "utl/Logging.h"
#include <cstdlib>

namespace dhsm
{

class HardwareFailure : public std::exception
{
public:
	HardwareFailure() {}
	~HardwareFailure() throw() {}
};

struct CruiseControlData : StateData
{
	int speed;

	CruiseControlData() : speed(0) {}

	void onStart()
	{
		LOG_INFO() << "starting cruise control";
	}
	void setSpeed(int newSpeed)
	{
		speed = newSpeed;
		LOG_INFO() << "setting speed at " << speed << " km/h";
	}
	void resumeSpeed()
	{
		LOG_INFO() << "resuming at speed " << speed << " km/h";
	}
};

struct CruiseData : StateData
{
	CruiseData()
	{
		LOG_INFO() << "building cruise data";
	}
	~CruiseData()
	{
		LOG_INFO() << "destroying cruise data";
	}
	void onCruiseBegin()
	{
		LOG_INFO() << "speed locked at " << dynamic_cast<CruiseControlData*>(context_)->speed << " km/h";
	}
};

// =============================================================================
// dhsm::CruiseControl class methods body
// -----------------------------------------------------------------------------

const char CruiseControl::EV_EN_START[] = "EV_EN_START";
const char CruiseControl::EV_EN_STOP[] = "EV_EN_STOP";
const char CruiseControl::ST_EN_OFF[] = "ST_EN_OFF";
const char CruiseControl::ST_EN_ON[] = "ST_EN_ON";

const char CruiseControl::EV_CC_START[] = "EV_CC_START";
const char CruiseControl::EV_CC_STOP[] = "EV_CC_STOP";
const char CruiseControl::ST_CC_OFF[] = "ST_CC_OFF";
const char CruiseControl::ST_CC_ON[] = "ST_CC_ON";
const char CruiseControl::ST_CC_DISABLED[] = "ST_CC_DISABLED";
const char CruiseControl::ST_CC_ENABLED[] = "ST_CC_ENABLED";

const char CruiseControl::EV_SET_SPEED[] = "EV_SET_SPEED";
const char CruiseControl::EV_RESUME_CC[] = "EV_RESUME_CC";

const char CruiseControl::EV_ACCELERATE[] = "EV_ACCELERATE";
const char CruiseControl::EV_BRAKE[] = "EV_BRAKE";

CruiseControl::CruiseControl()
	: hwFailure(false)
{
	StateMachinePtr sm(new StateMachine("SM_ENGINE"));
	sm->initial(ST_EN_OFF).actions().add(*this, &CruiseControl::onEngineInit);
	sm->transition(EV_EN_START, ST_EN_OFF, ST_EN_ON);
	sm->transition(EV_EN_STOP, ST_EN_ON, ST_EN_OFF);
	sm->findState(ST_EN_OFF).enterActions()
		.add(*this, &CruiseControl::onEngineStopped);
	sm->findState(ST_EN_ON).enterActions()
		.add(*this, &CruiseControl::onEngineStarted);
	sm->resetActions().add(*this, &CruiseControl::onEngineReset);

	StateMachine& engineOnSm = sm->state(ST_EN_ON).stateMachine("SM_CCONTROL");
	engineOnSm.state<CruiseControlData>(ST_CC_ON).enterActions()
		.add(&CruiseControlData::onStart);
	engineOnSm.initial(ST_CC_OFF).actions()
		.add(*this, &CruiseControl::onCruiseControlInit);
	engineOnSm.transition(EV_CC_START, ST_CC_OFF, ST_CC_ON);
	engineOnSm.transition(EV_CC_STOP, ST_CC_ON, ST_CC_OFF);
	engineOnSm.findState(ST_CC_OFF).enterActions()
		.add(*this, &CruiseControl::onCruiseStopped);
	engineOnSm.findState(ST_CC_ON).enterActions()
		.add(*this, &CruiseControl::onCruiseStarted);
	engineOnSm.resetActions()
		.add(*this, &CruiseControl::onCruiseControlReset);
	engineOnSm.errState().enterActions()
		.add(*this, &CruiseControl::onCruiseControlError);

	StateMachine& cruiseOnSm = engineOnSm.state(ST_CC_ON).stateMachine("SM_CRUISE");
	cruiseOnSm.state<CruiseData>(ST_CC_ENABLED);
	cruiseOnSm.initial(ST_CC_DISABLED).actions()
		.add(*this, &CruiseControl::onActiveCruiseInit);
	cruiseOnSm.transition(EV_SET_SPEED, ST_CC_DISABLED, ST_CC_ENABLED).actions()
		.add(&CruiseControlData::setSpeed)
		.add(*this, &CruiseControl::onSetSpeed);
	cruiseOnSm.transition(EV_RESUME_CC, ST_CC_DISABLED, ST_CC_ENABLED).actions()
		.add(&CruiseControlData::resumeSpeed);
	cruiseOnSm.transition(EV_ACCELERATE, ST_CC_ENABLED, ST_CC_DISABLED);
	cruiseOnSm.transition(EV_BRAKE, ST_CC_ENABLED, ST_CC_DISABLED);
	cruiseOnSm.findState(ST_CC_DISABLED).enterActions()
		.add(*this, &CruiseControl::onCruiseDisabled);
	cruiseOnSm.findState(ST_CC_ENABLED).enterActions()
		.add(*this, &CruiseControl::onCruiseEnabled)
		.add(&CruiseData::onCruiseBegin);
	cruiseOnSm.resetActions()
		.add(*this, &CruiseControl::onActiveCruiseReset);
	cruiseOnSm.errState().enterActions()
		.add(*this, &CruiseControl::onCruiseError);

	engine.setStateMachine(sm);
}

void CruiseControl::start()
{
	engine.start();
}

void CruiseControl::startEngine()
{
	LOG_INFO() << ">>>> start engine";
	engine.step(EV_EN_START);
}

void CruiseControl::stopEngine()
{
	LOG_INFO() << ">>>> stop engine";
	engine.step(EV_EN_STOP);
}

void CruiseControl::startCruise()
{
	LOG_INFO() << ">>>> start cruise";
	engine.step(EV_CC_START);
}

void CruiseControl::stopCruise()
{
	LOG_INFO() << ">>>> stop cruise";
	engine.step(EV_CC_STOP);
}

void CruiseControl::setSpeed(int speed)
{
	LOG_INFO() << ">>>> set speed";
	Message<int> msg(EV_SET_SPEED, speed);
	engine.step(msg);
}

void CruiseControl::resumeCruise()
{
	LOG_INFO() << ">>>> resume cruise";
	engine.step(EV_RESUME_CC);
}

void CruiseControl::accelerate()
{
	LOG_INFO() << ">>>> accelerate";
	engine.step(EV_ACCELERATE);
}

void CruiseControl::brake()
{
	LOG_INFO() << ">>>> brake";
	engine.step(EV_BRAKE);
}

void CruiseControl::setHwFailure(bool failure)
{
	hwFailure = failure;
}

bool CruiseControl::hasFailed()
{
	return engine.inErrorState();
}

void CruiseControl::recoverFailure()
{
	setHwFailure(false);
	engine.reinitiate(engine.currentStateMachine());
}

void CruiseControl::onEngineStopped()
{
	LOG_INFO() << "engine turned off";
}

void CruiseControl::onEngineStarted()
{
	LOG_INFO() << "engine turned on";
}

void CruiseControl::onCruiseStopped()
{
	LOG_INFO() << "cruise control turned off";
}

void CruiseControl::onCruiseStarted()
{
	LOG_INFO() << "cruise control turned on";
}

void CruiseControl::onCruiseEnabled()
{
	LOG_INFO() << "cruise control enabled";
}

void CruiseControl::onCruiseDisabled()
{
	LOG_INFO() << "cruise control disabled";
}

void CruiseControl::onSetSpeed(int speed)
{
	LOG_INFO() << "speed set at " << speed << " km/h";
}

void CruiseControl::onEngineInit()
{
	LOG_INFO() << "engine initiated";
}

void CruiseControl::onCruiseControlInit()
{
	if(hwFailure)
	{
		throw HardwareFailure();
	}
	LOG_INFO() << "cruise control initiated";
}

void CruiseControl::onActiveCruiseInit()
{
	if(hwFailure)
	{
		throw HardwareFailure();
	}
	LOG_INFO() << "new cruise initiated";
}

void CruiseControl::onEngineReset()
{
	LOG_INFO() << "engine reset!";
}

void CruiseControl::onCruiseControlReset()
{
	LOG_INFO() << "cruise control reset!";
}

void CruiseControl::onActiveCruiseReset()
{
	LOG_INFO() << "active cruise reset!";
}

void CruiseControl::onCruiseError()
{
	LOG_INFO() << "cruise ERROR";
	try
	{
		throw;
	}
	catch(const HardwareFailure& e)
	{
		throw;
	}
	catch(...)
	{
	}
}

void CruiseControl::onCruiseControlError()
{
	LOG_INFO() << "cruise control ERROR";
}

} // namespace dhsm
