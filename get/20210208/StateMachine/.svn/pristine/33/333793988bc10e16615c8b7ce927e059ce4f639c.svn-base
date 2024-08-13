// =============================================================================
// State.h
// -----------------------------------------------------------------------------
// class dhsm::State
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_State_h_INCLUDED
#define dhsm_State_h_INCLUDED

#include "Action.hpp"
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <string>
#include <memory>

namespace dhsm
{
// Forward declarations
class StateMachine;
class Transition;
typedef std::auto_ptr<StateMachine> StateMachinePtr;

/**
 * This class stores the static characteristics of a state within a state machine.
 *
 * It has an associated label for human intelligibility and display.
 * It has two list of actions set by the user to implement any specific behavior
 * when a transition causes the state machine to enter or exit this state.
 */
class State
{
public:
	/// This container class maps Events on Transitions for this State.
	typedef boost::ptr_map<Event, Transition> TransitionMap;

public:
	State(StateMachine* parentSm, const char* name);
	//@{
	const std::string& label() const { return label_; }
	void setLabel(const std::string& newLabel) { label_ = newLabel; }
	std::string path() const;
	void path(std::ostringstream& s) const;
	//@}
	template<class StateDataType>
	void setStateDataType();

	//@{
	Transition& addTransition(const Event& triggerEvent, State& targetState);
	Transition* getTransition(const Event& event);
	const TransitionMap& transitionMap() const { return transitionMap_; }
	ActionList& enterActions() { return enterActions_; }
	ActionList& exitActions() { return exitActions_; }
	bool isSensitive(const Event& event) const;
	//@}
	//@{
	StateMachine* parentStateMachine() { return parentStateMachine_; }
	const StateMachine* parentStateMachine() const { return parentStateMachine_; }
	bool isComposite() const { return stateMachine_.get() != 0; }
	StateMachine& stateMachine(const char* name = 0);
	StateMachine* getStateMachine() const { return stateMachine_.get(); }
	StateMachine& setStateMachine(std::auto_ptr<StateMachine> subSm);
	//@}

	StateDataPtr createStateData(StateData* context, Engine* engine);

protected:
	/// This is the type for methods responsible for state data instanciation
	typedef StateDataPtr (*FactoryMethod)(const StateData& base);

	template<class T>
	static StateDataPtr factoryMethod(const typename boost::enable_if<boost::is_convertible<StateData, const T&>, StateData>::type& base);

	template<class T>
	static StateDataPtr factoryMethod(const typename boost::disable_if<boost::is_convertible<StateData, const T&>, StateData>::type& base);

	std::string label_; ///< the name of the state.
	/// All transitions associated with this state. The map's key is the Event
	/// that triggers the Transition in this State.
	TransitionMap transitionMap_;
	/// The list of actions that should be executed when entering this state.
	ActionList enterActions_;
	/// The list of actions that should be executed when leaving this state.
	ActionList exitActions_;
	/// The nested state machine, if any.
	std::auto_ptr<StateMachine> stateMachine_;
	/// The state machine this state belongs to.
	StateMachine* parentStateMachine_;
	/// The factory method that can create the subclass of StateData associated
	// with this state.
	FactoryMethod createStateData_;
};


// =============================================================================
// Inline methods of class State
// -----------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const State& state)
{
	stream << state.label();
	return stream;
}

/**
 * Creates the StateData object associated with this state.
 */
inline StateDataPtr State::createStateData(StateData* context, Engine* engine)
{
	StateDataPtr sdata;
	if(createStateData_ != 0)
	{
		StateData base(this, context, engine);
		sdata.reset(createStateData_(base).release());
	}
	return sdata;
}

/**
 * Defines the class of the object whose lifetime should be bound to the
 * activation if this state.
 * The parameter class StateDataType must inherit from StateData.
 * Use the type "void" to specify that no object should be created (this is the
 * default).
 */
template<class StateDataType>
void State::setStateDataType()
{
	createStateData_ = &State::factoryMethod<StateDataType>;
}

/**
 * Implements the instanciation of the object bound to this state.
 */
template<class T>
StateDataPtr State::factoryMethod(const typename boost::enable_if<boost::is_convertible<StateData, const T&>, StateData>::type& base)
{
	return StateDataPtr(new T(base));
}

template<class T>
StateDataPtr State::factoryMethod(const typename boost::disable_if<boost::is_convertible<StateData, const T&>, StateData>::type& base)
{
	StateDataPtr data(new T);
	data->setState(base.state());
	data->setContext(base.context());
	data->setEngine(base.engine());
	return data;
}

} // namespace dhsm

#endif // dhsm_State_h_INCLUDED
