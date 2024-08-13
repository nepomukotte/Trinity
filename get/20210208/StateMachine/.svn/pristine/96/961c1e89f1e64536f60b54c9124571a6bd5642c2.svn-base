// =============================================================================
// RunControlStateMachine.h
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

#ifndef dhsm_RunControlStateMachine_h_INCLUDED
#define dhsm_RunControlStateMachine_h_INCLUDED

#include "dhsm/StateMachine.h"
#include "RunControlBehaviour.h"

namespace dhsm
{
/**
 * This class implements the execution logic of the run control test state-machine.
 */
class RunControlStateMachine : public StateMachine
{
public:
	class MyCriticalError;

	RunControlStateMachine(Behaviour behaviour);
	~RunControlStateMachine();

private:
	void onTurnedOn();
	void onTurnedOff();
	void onIdling();
	void onActivated();
	void onRunning();
	void onPaused();
	void onError();

	Behaviour behaviour;
};

class RunControlStateMachine::MyCriticalError : public CriticalError
{
public:
	MyCriticalError();
	~MyCriticalError() throw();

	void writeMessage(std::ostream& stream) const;
};

} // namespace dhsm

#endif // dhsm_RunControlStateMachine_h_INCLUDED

