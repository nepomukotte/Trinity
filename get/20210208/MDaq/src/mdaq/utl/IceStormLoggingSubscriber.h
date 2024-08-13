/**
 * @file IceStormLoggingSubscriber.h
 * @date 6 déc. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the MDaq software project.
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

#ifndef mdaq_utl_IceStormLoggingSubscriber_h_INCLUDED
#define mdaq_utl_IceStormLoggingSubscriber_h_INCLUDED

#include "mdaq/utl/LoggingService.h"
#include "mdaq/DefaultPortNums.h"
#include "utl/net/SocketAddress.h"
#include ICE_HEADER(Communicator.h)
#include <IceStorm/IceStorm.h>

namespace mdaq
{
namespace utl
{
/*
 * An implementation of the IceStorm remote logging contract.
 */
class LoggingServiceI : virtual public LoggingService
{
public:
	virtual void log(::mdaq::utl::LogSeverity level, const ::std::string & source,
			const ::std::string & msg, const ::Ice::Current & = ::Ice::Current());
static IceStorm::TopicPrx getTopic(const std::string & topicName = "Logger");
};

/**
 * A class to subscribe to log messages sent via IceStorm.
 */
class IceStormLoggingSubscriber
{
public:
	IceStormLoggingSubscriber(const ::utl::net::SocketAddress & address, const std::string & topicName = "Logger");
	IceStormLoggingSubscriber(const ::utl::net::IpAddress & address, const std::string & topicName = "Logger");
	IceStormLoggingSubscriber(const uint16_t port, const std::string & topicName = "Logger");
	IceStormLoggingSubscriber(const std::string & topicName = "Logger");
	~IceStormLoggingSubscriber();
	bool isValid() const { return isValid_; }
private:
	Ice::ObjectPrx createSubscriber(uint16_t port);
	Ice::ObjectPrx createSubscriber(const ::utl::net::SocketAddress & addr);
private:
	IceStorm::TopicPrx topic;
	Ice::ObjectPrx subscriber;
	bool isValid_;
};

} /* namespace utl */
} /* namespace mdaq */
#endif /* mdaq_utl_IceStormLoggingSubscriber_h_INCLUDED */
