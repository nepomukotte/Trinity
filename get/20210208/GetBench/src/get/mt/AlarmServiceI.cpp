/**
 * @file AlarmServiceI.cpp
 * @date Jan 14, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: AlarmServiceI.cpp 1244 2014-04-17 09:20:40Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the GetBench software project.
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

#include "AlarmServiceI.h"
#include "get/cobo/CtrlNodeI.h"
#include "utl/Logging.h"
#include "utl/BitFieldHelper.hpp"
#include "boost/cstdint.hpp"

namespace get
{
namespace mt
{

AlarmServiceI::AlarmServiceI(cobo::CtrlNodeI & ctrlNode) : cobo_(ctrlNode)
{
}

void AlarmServiceI::subscribe(const ::std::string & id, const ::get::mt::AlarmCallbackPrx & cb, const ::Ice::Current &)
{
	IceUtil::Mutex::Lock lock(mutex_);
	if (cb)
	{
		LOG_DEBUG() << "Client '" << id << "' subscribed to alarm service";
		subscribers_[id] = cb;
	}
}

void AlarmServiceI::sendAsadAlarm(const ::std::string & coboId, ::Ice::Long alarmBits, const ::Ice::Current &)
{
	IceUtil::Mutex::Lock lock(mutex_);
	//LOG_DEBUG() << "Sending AsAd alarm to all subscribers...";
	AlarmCallbackMap::iterator s = subscribers_.begin();
	while (s != subscribers_.end())
	{
		try
		{
			s->second->sendAsadAlarm(coboId, alarmBits);
			++s;
		}
		catch (const Ice::LocalException & e)
		{
			LOG_ERROR() << "Could not send AsAd alarm to client '" << s->first << "': " << e.ice_name();
			subscribers_.erase(s++);
		}
	}
	// Pulse CoBo's red 'B' LED's on
	const uint8_t asadAlarmMask = static_cast<uint8_t>(alarmBits);
	for (size_t asadId=0; asadId < 4u; ++asadId)
	{
		if (::utl::BitFieldHelper< uint8_t >::getBit(asadAlarmMask, asadId))
			cobo_.modifyLED(cobo::LedB, asadId, cobo::LedFastPulse);
	}
}

void AlarmServiceI::unsubscribe(const std::string & id, const Ice::Current &)
{
	IceUtil::Mutex::Lock lock(mutex_);
	subscribers_.erase(id);
}

void AlarmServiceI::reset(const ::Ice::Current &)
{
	IceUtil::Mutex::Lock lock(mutex_);
	subscribers_.clear();
}

} /* namespace mt */
} /* namespace get */
