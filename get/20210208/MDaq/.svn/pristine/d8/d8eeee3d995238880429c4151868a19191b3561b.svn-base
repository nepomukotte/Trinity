/**
 * @file Log4cxxLoggingBackend.h
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

#ifndef mdaq_utl_Log4cxxLoggingBackend_h_INCLUDED
#define mdaq_utl_Log4cxxLoggingBackend_h_INCLUDED

#include "utl/Logging.h"
#include <log4cxx/logger.h>
#include <sstream>

namespace mdaq
{
namespace utl
{
/**
 * This is a logging backend that uses the Apache log4cxx library.
 */
class Log4cxxLoggingBackend : public ::utl::AbstractLoggingBackend
{
public:
	Log4cxxLoggingBackend(const std::string & name="mdaq");
	virtual ~Log4cxxLoggingBackend();
	std::ostream& beginLog(::utl::SeverityLevel level, const char* source, unsigned line);
	void endLog();
private:
	std::string loggerName(const std::string & source) const;
	std::ostringstream stream;
	::utl::SeverityLevel currentLevel;
	std::string baseName; ///< Logger base name.
	::log4cxx::LoggerPtr logger;
};

} /* namespace utl */
} /* namespace mdaq */
#endif /* mdaq_utl_Log4cxxLoggingBackend_h_INCLUDED */
