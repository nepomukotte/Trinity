// =============================================================================
// Engine.h
// -----------------------------------------------------------------------------
// class dhsm::Engine
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_Engine_h_INCLUDED
#define dhsm_Engine_h_INCLUDED

#include "StateMachine.h"
#include "Event.h"
#include <set>
#include <string>
#include <climits>

namespace dhsm
{
/**
 * This class lays the bases of the dynamic behaviour of a state machine.
 * It only supports direct event processing (no queue).
 *
 * The engine uses the static description of a state machine (ie: a top-level
 * StateMachine object) and "gives it life".
 *
 * First the user should build a state machine, then associate it with an engine
 * using the setStateMachine method. Next, he should call initiate() to trigger
 * the first transition, and finally send events with the step() method to make
 * it run.
 *
 * @author fchateau
 */
class Engine
{
public:
	// enums
	enum ProcessingResult { Processed = 0, Ignored };

	Engine();
	virtual ~Engine();
	//@{
	StateMachine* stateMachine() { return stateMachine_.get(); }
	virtual void setStateMachine(StateMachinePtr sm);
	StateMachinePtr releaseStateMachine();
	//@}
	//@{
	void setLabel(const std::string& newLabel) { label_ = newLabel; }
	const std::string& label() const { return label_; }
	//@{
	void setTopContext(StateDataPtr ctx);
	StateData* getTopContext() { return topContext_.get(); }
	template<class T> T& topContext() { return dynamic_cast<T&>(*topContext_); }
	//@}
	//@{
	State& currentState(StateMachine& sm) { return *sm.currentState_; }
	State& currentState(int level = INT_MAX);
	StateMachine& currentStateMachine(int level = INT_MAX) { return *currentState(level).parentStateMachine(); }
	StateMachine* searchCurrentStateMachine(const std::string& smLabel);
	StateMachine& findCurrentStateMachine(const std::string& smLabel);
	bool isInitiated(StateMachine& sm) const  { return sm.currentState_ != sm.OFF_STATE; }
	bool inErrorState(StateMachine& sm) const { return sm.currentState_ == sm.ERR_STATE; }
	bool inErrorState();
	//@}
	//@{
	ProcessingResult start();
	ProcessingResult halt(StateMachine& sm);
	ProcessingResult reinitiate(StateMachine& sm);
	ProcessingResult step(const Event& event);
	//@}
protected:
	//@{
	virtual ProcessingResult doStart();
	virtual ProcessingResult doHalt(StateMachine& sm);
	virtual ProcessingResult doReinitiate(StateMachine& sm);
	virtual ProcessingResult crossTransition(StateMachine& sm, const Event& event, Transition& transition);
	//@}
	//@{
	static void clearStateData(StateMachine& sm);
	void createStateData(StateMachine& sm);
	//@}
	class StateDataDisposer;

protected:
	StateMachinePtr stateMachine_; ///< the top level state machine controlled by this engine.
	StateDataPtr topContext_;
	std::string label_;

private:
	//@{
	ProcessingResult initiate(StateMachine& sm);
	ProcessingResult process(StateMachine& sm, const Event& event);
	ProcessingResult reset(StateMachine& sm);
	//@}
	//@{
	ProcessingResult onCurrentStateEnter(StateMachine& sm, const Event& event);
	void onCurrentStateExit(StateMachine& sm);
	void onCurrentStateFailsafeExit(StateMachine& sm);
	void onTransitionTrigger(Transition& transition, StateData* closureData, const Event& event);
	void onStateMachineReset(StateMachine& sm, StateData* closureData);
	//@}

	inline void initSupportedEvents();
	inline void checkSupported(const Event& event);
#ifndef NDEBUG
	void listSupportedEvents(StateMachine& sm);

	typedef std::set<Event> EventSet;
	EventSet supportedEvents;
#endif
};

typedef std::auto_ptr<Engine> EnginePtr;

}

#endif // dhsm_Engine_h_INCLUDED
