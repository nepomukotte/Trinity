/**
 * @file Log4cxxLoggingBackend.cpp
 * @date 11 déc. 2012
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

#include "Log4cxxLoggingBackend.h"

namespace mdaq
{
namespace utl
{
//__________________________________________________________________________________________________
/**
 * Constructs a logging backend that transfers messages to an Apache log4cxx logger.
 * @param name Name of thr root logger.
 */
Log4cxxLoggingBackend::Log4cxxLoggingBackend(const std::string & name) :
		currentLevel(::utl::info), baseName(name), logger(log4cxx::Logger::getLogger(name))
{
}
//__________________________________________________________________________________________________
Log4cxxLoggingBackend::~Log4cxxLoggingBackend()
{
}
//__________________________________________________________________________________________________
std::string Log4cxxLoggingBackend::loggerName(const std::string & source) const
{
	return baseName + '.' + ::utl::xtract_base(source);
}
//__________________________________________________________________________________________________
std::ostream & Log4cxxLoggingBackend::beginLog(::utl::SeverityLevel level, const char* source, unsigned)
{
	logger = log4cxx::Logger::getLogger(loggerName(source));
	currentLevel = level;
	stream.str(std::string());
	return stream;
}
//__________________________________________________________________________________________________
void Log4cxxLoggingBackend::endLog()
{
	switch(currentLevel)
	{
	case ::utl::fatal:
		LOG4CXX_FATAL(logger, stream.str());
		break;
	case ::utl::error:
		LOG4CXX_ERROR(logger, stream.str());
		break;
	case ::utl::warn:
		LOG4CXX_WARN(logger, stream.str());
		break;
	case ::utl::debug:
		LOG4CXX_DEBUG(logger, stream.str());
		break;
	case ::utl::info:
	default:
		LOG4CXX_INFO(logger, stream.str());
		break;
	}
}
//__________________________________________________________________________________________________
} /* namespace utl */
} /* namespace mdaq */
