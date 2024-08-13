// =============================================================================
// StateMachine.h
// -----------------------------------------------------------------------------
// class dhsm::StateMachine
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_StateMachine_h_INCLUDED
#define dhsm_StateMachine_h_INCLUDED

#include "Transition.h"
#include "State.h"
#include "Event.h"
#include "Action.hpp"
#include "Exception.h"
#include <map>
#include <ostream>
#include <string>
#include <memory>

namespace dhsm
{
/**
 * This class represents a state-machine like those in UML, and provides methods
 * to set it up at run-time.
 *
 * A state machine is created with its initial pseudostate and final pseudostate
 * regrouped in a single state called the OFF State (OFF_STATE).
 * The OFF state corresponds to an inactive state machine (not yet initiated):
 * in this state, no event can be processed by the engine (Engine and its
 * subclasses).
 * A state machine can only be initiated if it has an initial transition. It
 * should be defined by calling the initial() method.
 *
 * It is possible to nest state machines in states in order to build a
 * hierarchical state machine.
 *
 * This class primarily focuses on describing the static structure of a state
 * machine. Most of the runtime support is provided by the engine classes.
 * Nevertheless, it is not a pure static representation, as it embeds some state
 * variables about the statemachine current state.
 * Although in theory, it could have been moved to the engine classes, in
 * practice, it would have introduced a non-trivial complexity and runtime
 * overhead (in particular to implement the concept of history), and the
 * benefits were too few.
 *
 * @author fchateau, sanvar
 */
class StateMachine
{
public:
	// Friend classes
	friend class Engine; ///< allows Engine to access currentState_ and currentStateData_
	friend class State; /// < allows State to modify the parentState

	// Typedefs
	typedef std::map<std::string, State*> StateMap;

	// Enums
	enum HistoryMode { HISTORY, NO_HISTORY };

	// Constants
	static const char OFF_STATE_LABEL[];
	static const char ERR_STATE_LABEL[];
	static const char INITIATE_EVENT_LABEL[];

public:
	StateMachine(const char* label = "STATEMACHINE", State* parentState = 0);
	virtual ~StateMachine();
	//@{
	std::string path() const;
	const std::string& label() const { return label_; }
	void setLabel(const std::string& newLabel) { label_ = newLabel; }
	HistoryMode historyMode() const { return historyMode_; }
	ActionList& resetActions() { return resetActions_; }
	State* parentState() { return parentState_; }
	const State* parentState() const { return parentState_; }
	StateMachine* parentStateMachine() { return parentState_ != 0 ? parentState_->parentStateMachine() : 0; }
	StateMachine& topParentStateMachine();
	StateData* parentCurrentStateData();
	StateData* currentStateData() { return currentStateData_.get(); }
	int nestingLevel() const;
	//@}
	//@{
	template<class T>
	State& state(const char* stateLabel);
	State& state(const char* stateLabel);
	State& findState(const std::string& stateLabel);
	State* searchState(const std::string& stateLabel);
	State& offState() { return *OFF_STATE; }
	State& errState() { return *ERR_STATE; }
	//@}
	//@{
	Transition& transition(
		const char* event,
		const char* sourceState,
		const char* targetState
	);
	Transition& transition(
		const Event& event,
		State& sourceState,
		State& targetState
	);
	Transition& initial(const char* entryState, HistoryMode history = NO_HISTORY);
	Transition& initial(State& entryState, HistoryMode history = NO_HISTORY);
	//@}

	#ifdef UTL_LOGGING_ENABLED
	void dump() const;
	#endif

private:
	std::string label_;

	State* OFF_STATE; ///< The state machine's OFF state.
	State* ERR_STATE; ///< The state machine's error state.
	HistoryMode historyMode_; ///< Indicates whether the state of the state machine is kept when leaving the composite state.

	State* parentState_; ///< The parent composite state. Null if this is the top state machine.
	StateMap stateMap_; ///< Maps State labels to States.
	ActionList resetActions_; ///< The list of actions to execute when the state machine gets reset.

	State* currentState_; ///< The state machine's current state.
	StateDataPtr currentStateData_;
	bool currentStateIsReady_;
};

// =============================================================================
// dhsm::StateMachine class inline methods body
// -----------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const StateMachine& sm)
{
	stream << '#' << sm.label();
	return stream;
}

/**
 * Adds a new State to the StateMachine's set of States.
 * It also defines the class for the associated data object that gets
 * instanciated when the state becomes active and destroyed when it becomes
 * inactive.
 * If a state with the same label already exists, the data type associate with
 * the state will be modified.
 */
template<class T>
State& StateMachine::state(const char* stateLabel)
{
	State& st = state(stateLabel);
	st.setStateDataType<T>();
	return st;
}

} // namespace dhsm

#endif // dhsm_Engine_h_INCLUDED
