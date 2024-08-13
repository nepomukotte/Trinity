// =============================================================================
// StateData.h
// -----------------------------------------------------------------------------
// class dhsm::StateData
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_StateData_h_INCLUDED
#define dhsm_StateData_h_INCLUDED

#include <memory>

namespace dhsm
{
// forward declarations
class State;
class Engine;

/**
 * This is the base class for user-defined classes associated to a state and
 * holding application specific data whose lifetime is bound to the state activation.
 */
class StateData
{
public:
	StateData() : state_(0), context_(0), engine_(0) {}
	StateData(State* state, StateData* context, Engine* engine) : state_(state), context_(context), engine_(engine) {}
	virtual ~StateData() {}

	State* state() const { return state_; }
	StateData* context() const { return context_; }
	Engine* engine() const { return engine_; }

	void setState(State* st) { state_ = st; }
	void setContext(StateData* ctxt) { context_ = ctxt; }
	void setEngine(Engine* engn) { engine_ = engn; }

protected:
	State* state_; ///< the state this object is associated to.
	StateData* context_;  ///< the state data associated with the parent state.
	Engine* engine_; ///< the engine that runs the state machine and that created this object.
};

typedef std::auto_ptr<StateData> StateDataPtr;

}

#endif // dhsm_StateData_h_INCLUDED
