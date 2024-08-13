/**
 * @file IceStormLoggingBackend.cpp
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

#include "IceStormLoggingBackend.h"
#include "IceStormLoggingSubscriber.h"
#include "utl/Logging.h"

namespace mdaq
{
namespace utl
{
//__________________________________________________________________________________________________
LogSeverity IceStormLoggingBackend::loc2iceSeverityMap[] = {debug, info, warn, error, fatal};
//__________________________________________________________________________________________________
IceStormLoggingBackend::IceStormLoggingBackend(const std::string & topicName) : topic(0), loggingService(0)
{
	topic = LoggingServiceI::getTopic(topicName);
	createPublisher();
}
//__________________________________________________________________________________________________
void IceStormLoggingBackend::createPublisher()
{
	if (topic)
	{
		Ice::ObjectPrx pub = topic->getPublisher()->ice_oneway();
		loggingService = LoggingServicePrx::uncheckedCast(pub);
		LOG_DEBUG() << "Successfully connected to IceStorm logging topic.";
	}
}
//__________________________________________________________________________________________________
IceStormLoggingBackend::~IceStormLoggingBackend()
{
}
//__________________________________________________________________________________________________
std::ostream & IceStormLoggingBackend::beginLog(::utl::SeverityLevel level, const char* source, unsigned line)
{
	if (! loggingService)
	{
		return fallback.beginLog(level, source, line);
	}
	else
	{
		currentLevel = level;
		currentSource = source;
		currentLine = line;
		stream.str(std::string());
		return stream;
	}
}
//__________________________________________________________________________________________________
void IceStormLoggingBackend::endLog()
{
	if (! loggingService)
	{
		fallback.endLog();
	}
	else
	{
		try
		{
			// Local logger
			fallback.beginLog(currentLevel, currentSource, currentLine) << stream.str();
			fallback.endLog();

			// Remote loggers
			loggingService->log(loc2iceSeverityMap[currentLevel], ::utl::xtract_base(currentSource), stream.str());
		}
		catch (const Ice::Exception& e)
		{
			loggingService = 0;
			LOG_WARN() << "Lost connection to remote logging service. Falling back on local console logger: " << e.what();
		}
	}
}
//__________________________________________________________________________________________________
} /* namespace utl */
} /* namespace mdaq */
