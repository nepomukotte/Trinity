// =============================================================================
// CruiseControl.h
// -----------------------------------------------------------------------------
// class dhsm::CruiseControl
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_CruiseControl_h_INCLUDED
#define dhsm_CruiseControl_h_INCLUDED

#include "dhsm/Engine.h"

namespace dhsm
{
class CruiseControl
{
public:
	CruiseControl();

	void start();
	void startEngine();
	void stopEngine();
	void startCruise();
	void stopCruise();
	void setSpeed(int speed);
	void resumeCruise();
	void accelerate();
	void brake();

	void setHwFailure(bool failure);
	bool hasFailed();
	void recoverFailure();

private:
	static const char EV_EN_START[];
	static const char EV_EN_STOP[];
	static const char ST_EN_OFF[];
	static const char ST_EN_ON[];

	static const char EV_CC_START[];
	static const char EV_CC_STOP[];
	static const char ST_CC_OFF[];
	static const char ST_CC_ON[];
	static const char ST_CC_DISABLED[];
	static const char ST_CC_ENABLED[];

	static const char EV_SET_SPEED[];
	static const char EV_RESUME_CC[];

	static const char EV_ACCELERATE[];
	static const char EV_BRAKE[];

	void onEngineStopped();
	void onEngineStarted();
	void onCruiseStopped();
	void onCruiseStarted();
	void onCruiseEnabled();
	void onCruiseDisabled();

	void onSetSpeed(int speed);

	void onEngineInit();
	void onCruiseControlInit();
	void onActiveCruiseInit();

	void onEngineReset();
	void onCruiseControlReset();
	void onActiveCruiseReset();

	void onCruiseError();
	void onCruiseControlError();

	Engine engine;
	bool hwFailure;
};

} // namespace dhsm

#endif // dhsm_CruiseControl_h_INCLUDED
