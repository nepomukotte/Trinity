// =============================================================================
// Chrono.h
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

#ifndef dhsm_Chrono_h_INCLUDED
#define dhsm_Chrono_h_INCLUDED

#include "dhsm/Engine.h"
#include <string>

namespace dhsm
{
class Chrono : public Engine
{
public:
	Chrono();
	~Chrono();

	int getCount();
	void dump();

	static const char ST_IDLING[];
	static const char ST_ACTIVE[];
	static const char ST_PAUSED[];
	static const char EV_TOGGLE[];
	static const char EV_TRIGGR[];
	static const char EV_SUSPND[];
	static const char EV_RESUME[];
	static const char EV_DSABLE[];

private:
	void onStart();
	void onStop();
	void onPause();
	void onResume();
	void onLaunch();
	void onIncrCount();
	void onIncrMiss();

private:
	int count;
	int miss;
};

} // namespace dhsm


#endif // Chrono_h_INCLUDED
