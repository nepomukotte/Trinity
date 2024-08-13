/**
 * @file CompositeEngine.cpp
 * @date 27 mai 2013
 * @author fchateau
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Frederic Chateau
 *
 * This file is part of the StateMachine software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "CompositeEngine.h"

namespace dhsm
{

CompositeEngine::CompositeEngine()
{
}

CompositeEngine::~CompositeEngine()
{
}

void CompositeEngine::appendEngine(Engine* engine)
{
	engines_.push_back(engine);
}

void CompositeEngine::clearEngines()
{
	engines_.clear();
}

void CompositeEngine::setEventSequence(const std::string& event, bool sequenceForward)
{
	sequences_.insert(std::make_pair(event, sequenceForward));
}

CompositeEngine::ProcessingResult
CompositeEngine::doStart()
{
	for(EngineList::iterator engine = engines_.begin(), end = engines_.end(); engine != end; ++engine)
	{
		(*engine)->start();
	}
	return Engine::doStart();
}

CompositeEngine::ProcessingResult
CompositeEngine::doReinitiate(StateMachine& sm)
{
	for(EngineList::iterator engine = engines_.begin(), end = engines_.end(); engine != end; ++engine)
	{
		(*engine)->reinitiate((*engine)->findCurrentStateMachine(sm.label()));
	}
	return Engine::doReinitiate(sm);
}

CompositeEngine::ProcessingResult
CompositeEngine::doHalt(StateMachine& sm)
{
	for(EngineList::iterator engine = engines_.begin(), end = engines_.end(); engine != end; ++engine)
	{
		(*engine)->halt((*engine)->findCurrentStateMachine(sm.label()));
	}
	return Engine::doHalt(sm);
}

CompositeEngine::ProcessingResult
CompositeEngine::crossTransition(StateMachine& sm, const Event& event, Transition& transition)
{
	EventSequenceMap::const_iterator seq = sequences_.find(event.label());
	if(seq == sequences_.end())
	{
		throw StaticError(StaticError::UnsupportedEvent, event.label());
	}
	if(seq->second == true) // forward sequence
	{
		for(EngineList::iterator engine = engines_.begin(), end = engines_.end();
			engine != end; ++engine)
		{
			(*engine)->step(event);
		}
	}
	else // backward sequence
	{
		for(EngineList::reverse_iterator engine = engines_.rbegin(), end = engines_.rend();
			engine != end; ++engine)
		{
			(*engine)->step(event);
		}
	}
	return Engine::crossTransition(sm, event, transition);
}

} /* namespace dhsm */
