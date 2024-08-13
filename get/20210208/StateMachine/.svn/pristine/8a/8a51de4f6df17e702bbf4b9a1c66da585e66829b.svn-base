// =============================================================================
// EngineServant.h
// -----------------------------------------------------------------------------
// class dhsm::EngineServant
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_EngineServant_h_INCLUDED
#define dhsm_EngineServant_h_INCLUDED

#include "IEngineServant.h"
#include "Engine.h"

namespace dhsm
{
/**
 * This class is the implementation of the ICE servant for the state machine engine.
 * It awaits commands from a remote engine master and execute them locally on
 * its associated engine.
 */
class EngineServant : public IEngineServant
{
public:
	EngineServant();
	virtual ~EngineServant();

	void setEngine(EnginePtr& newEngine);
	Engine& engine() { return *engine_; }
	const Engine& engine() const { return *engine_; }

protected:
	//@{
	virtual int start(const Ice::Current&);
	virtual int reinitiate(const std::string& stateMachine, const Ice::Current&);
	virtual int halt(const std::string& stateMachine, const Ice::Current&);
	//@}
	//@{
	virtual int stepEvent(const std::string& event, const Ice::Current&);
	virtual int stepStringMsg(const std::string& event, const std::string& param, const Ice::Current&);
	virtual int stepIntMsg(const std::string& event, int param, const Ice::Current&);
	//@}
	virtual void handleExceptions();

private:
	EnginePtr engine_;	///< the local engine to which commands are forwarded
};

}

#endif // dhsm_EngineServant_h_INCLUDED
