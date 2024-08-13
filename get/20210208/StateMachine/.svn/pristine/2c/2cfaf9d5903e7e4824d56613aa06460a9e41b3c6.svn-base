// =============================================================================
// SmTest.cpp
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#define USE_LOG

#include "dhsm/Camera.h"
#include "dhsm/Chrono.h"
#include "dhsm/CruiseControl.h"
#include "utl/Logging.h"

using namespace dhsm;

void testCamera1()
{
	LOG_INFO() << "======== Camera State Machine Test";
	Camera camera;
	camera.start();
	camera.pressShutterHalf();
	camera.releaseShutter();
	camera.pressConfigure();
	camera.pressShutterHalf();
	camera.pressConfigure();
	camera.pressShutterHalf();
	camera.pressShutterFull();
	camera.releaseShutter();
	LOG_INFO() << ">>>> Camera test finished";
}

void testChrono1()
{
	LOG_INFO() << "======== Chrono State Machine Test";
	Chrono chrono;
	chrono.start();

	chrono.step(Chrono::EV_TOGGLE);
	chrono.step(Chrono::EV_SUSPND);
	chrono.step(Chrono::EV_RESUME);
	chrono.step(Chrono::EV_SUSPND);
	chrono.step(Chrono::EV_TOGGLE);
	chrono.step(Chrono::EV_DSABLE);
}

void testCruiseControl1()
{
	LOG_INFO() << "======== Cruise Control State Machine Test";
	CruiseControl c;
	c.start();
	c.startEngine();
	c.startCruise();
	c.setSpeed(70);
	c.accelerate();
	c.resumeCruise();
	c.brake();
	c.accelerate();
	c.resumeCruise();
	c.stopCruise();
	c.resumeCruise();
	c.brake();

	c.setHwFailure(true);
	bool failed = false;
	try
	{
		c.startCruise();
	}
	catch(const TransitionError& e)
	{
		c.recoverFailure();
		failed = true;
	}
	if(!failed)
	{
		LOG_INFO() << "!!!!! the state machine should have failed !!!!!";
	}

	c.startCruise();

	c.setSpeed(90);
	c.stopEngine();
	LOG_INFO() << ">>>> Cruise control test finished";
}

int main()
{
	try
	{
		testChrono1();
		testCruiseControl1();
		testCamera1();
		return 0;
	}
	catch(const std::exception& e)
	{
		LOG_ERROR() << "**** " << e.what();
	}
	catch(const char* e)
	{
		LOG_ERROR() << "**** " << e;
	}
	catch(...)
	{
		LOG_ERROR() << "**** " << "Unknown exception";
	}
	return -1;
}
