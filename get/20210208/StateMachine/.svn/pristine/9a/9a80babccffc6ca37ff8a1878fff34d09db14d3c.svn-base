// =============================================================================
// Camera.h
// -----------------------------------------------------------------------------
// class dhsm::Camera
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_Camera_h_INCLUDED
#define dhsm_Camera_h_INCLUDED

#include "dhsm/Engine.h"

namespace dhsm
{
class Camera
{
public:
	Camera();

	void start();
	void pressShutterHalf();
	void pressShutterFull();
	void releaseShutter();
	void pressConfigure();

private:
	Engine engine;
};

} // namespace dhsm

#endif // dhsm_Camera_h_INCLUDED
