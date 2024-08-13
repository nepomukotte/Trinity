/**
 * @file IceStormLoggingSubscriber.cpp
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

#include "IceStormLoggingSubscriber.h"
#include "utl/net/SocketAddress.h"
#include "utl/Logging.h"
#include "mdaq/utl/Server.h"

namespace mdaq
{
namespace utl
{
//__________________________________________________________________________________________________
/**
 * Transfers remote log messages to the local logging backend.
 */
void LoggingServiceI::log(::mdaq::utl::LogSeverity level, const ::std::string & source,
		const ::std::string & msg, const ::Ice::Current & )
{
	switch(level)
	{
	case fatal:
		LOG_FATAL() << '[' << source << "] " << msg;
		break;
	case error:
		LOG_ERROR() << '[' << source << "] " << msg;
		break;
	case warn:
		LOG_WARN() << '[' << source << "] " << msg;
		break;
	case info:
		LOG_INFO() << msg;
		break;
	case debug:
		LOG_DEBUG() << msg;
		break;
	default:
		break;
	}

}
//__________________________________________________________________________________________________
/**
 * Retrieves (or creates) the logging topic from the IceStorm topic manager.
 * @return Returns a proxy for the topic.
 */
IceStorm::TopicPrx LoggingServiceI::getTopic(const std::string & topicName)
{
	IceStorm::TopicPrx loggingTopic=0;
	Ice::CommunicatorPtr communicator = ::mdaq::utl::Server::ic();

	IceStorm::TopicManagerPrx topicManager = 0;
	try
	{
		Ice::ObjectPrx obj = communicator->propertyToProxy("TopicManager.Proxy");
		if (not obj)
		{
			LOG_WARN() << "Failed to retrieve logging topic '" << topicName
					<< "'. Could not find proxy for IceStorm topic manager.";
			return loggingTopic;
		}
		topicManager = IceStorm::TopicManagerPrx::checkedCast(obj);
	}
	catch (const Ice::LocalException & e)
	{
		LOG_WARN() << "Failed to retrieve logging topic '" << topicName
				<< "'. Could not connect to IceStorm topic manager (" << e.ice_name() << ").";
		return loggingTopic;
	}

    size_t numTrials = 2;
	while (!loggingTopic and numTrials > 0)
	{
		try
		{
			loggingTopic = topicManager->retrieve(topicName);
		}
		catch (const IceStorm::NoSuchTopic&)
		{
			try
			{
				loggingTopic = topicManager->create(topicName);
			}
			catch (const IceStorm::TopicExists&)
			{
				// Another client created the topic.
			}
		}
		numTrials--;
	}
	if (loggingTopic)
	{
		LOG_DEBUG() << "Successfully retrieved logging topic '" << topicName << "'.";
	}
	else
	{
		LOG_WARN() << "Failed to retrieve logging topic '" << topicName << "'.";
	}
	return loggingTopic;
}
//__________________________________________________________________________________________________
IceStormLoggingSubscriber::IceStormLoggingSubscriber(const ::utl::net::IpAddress & address, const std::string & topicName)
	: topic(0), subscriber(0), isValid_(false)
{
	// Retrieve topic
	topic = LoggingServiceI::getTopic(topicName);

	// Try 10 successive ports for creating subscriber
	if (topic)
	{
		uint16_t numTrials = 0;
		uint16_t port = mdaq::Default::remoteLoggerPortNum;
		while (!subscriber and numTrials < 10)
		{
			subscriber = createSubscriber(::utl::net::SocketAddress(address, port + numTrials));
			numTrials++;
		}
	}
	if (subscriber)
	{
		IceStorm::QoS qos;
		try
		{
			topic->subscribeAndGetPublisher(qos, subscriber);
			isValid_ = true;
		}
		catch (const IceStorm::AlreadySubscribed & e)
		{
			LOG_WARN() << "Error subscribing to topic '" << topicName << "': " << e.what();
		}
	}
}
//__________________________________________________________________________________________________
IceStormLoggingSubscriber::IceStormLoggingSubscriber(const ::utl::net::SocketAddress & address, const std::string & topicName)
	: topic(0), subscriber(0), isValid_(false)
{
	// Retrieve topic
	topic = LoggingServiceI::getTopic(topicName);

	// Try creating subscriber
	if (topic)
	{
		subscriber = createSubscriber(address);
	}
	if (subscriber)
	{
		IceStorm::QoS qos;
		try
		{
			topic->subscribeAndGetPublisher(qos, subscriber);
			isValid_ = true;
		}
		catch (const IceStorm::AlreadySubscribed & e)
		{
			LOG_WARN() << "Error subscribing to topic '" << topicName << "': " << e.what();
		}
	}
}
//__________________________________________________________________________________________________
IceStormLoggingSubscriber::IceStormLoggingSubscriber(const uint16_t port, const std::string & topicName)
	: topic(0), subscriber(0), isValid_(false)
{
	// Retrieve topic
	topic = LoggingServiceI::getTopic(topicName);

	// Try creating subscriber
	if (topic)
	{
		subscriber = createSubscriber(port);
	}
	if (subscriber)
	{
		IceStorm::QoS qos;
		try
		{
			topic->subscribeAndGetPublisher(qos, subscriber);
			isValid_ = true;
		}
		catch (const IceStorm::AlreadySubscribed & e)
		{
			LOG_WARN() << "Error subscribing to topic '" << topicName << "': " << e.what();
		}
	}
}
//__________________________________________________________________________________________________
IceStormLoggingSubscriber::IceStormLoggingSubscriber(const std::string & topicName)
	: topic(0), subscriber(0), isValid_(false)
{
	// Retrieve topic
	topic = LoggingServiceI::getTopic(topicName);

	// Try 10 successive ports for creating subscriber
	if (topic)
	{
		uint16_t numTrials = 0;
		uint16_t port = mdaq::Default::remoteLoggerPortNum;
		while (!subscriber and numTrials < 10)
		{
			subscriber = createSubscriber(port + numTrials);
			numTrials++;
		}
	}
	if (subscriber)
	{
		IceStorm::QoS qos;
		try
		{
			topic->subscribeAndGetPublisher(qos, subscriber);
			isValid_ = true;
		}
		catch (const IceStorm::AlreadySubscribed & e)
		{
			LOG_WARN() << "Error subscribing to topic '" << topicName << "': " << e.what();
		}
	}
}
//__________________________________________________________________________________________________
IceStormLoggingSubscriber::~IceStormLoggingSubscriber()
{
	if (topic and subscriber)
	{
		topic->unsubscribe(subscriber);
	}
}
//__________________________________________________________________________________________________
Ice::ObjectPrx IceStormLoggingSubscriber::createSubscriber(const ::utl::net::SocketAddress & addr)
{
	try
	{
		Server & server = Server::create(addr);
		std::string subscriberId = server.addServantWithUUID(new LoggingServiceI);
		Ice::ObjectPrx proxy = server.servantPrx(subscriberId)->ice_oneway();
		server.start();
		return proxy;
	}
	catch (const Ice::SocketException & e)
	{
		LOG_WARN() << "Could no create IceStorm logging subscriber using address " << addr << " ("
				<< e.ice_name() << ')';
	}
	return 0;
}
//__________________________________________________________________________________________________
Ice::ObjectPrx IceStormLoggingSubscriber::createSubscriber(uint16_t port)
{
	return createSubscriber(::utl::net::SocketAddress(0, port));
}
//__________________________________________________________________________________________________
} /* namespace utl */
} /* namespace mdaq */
