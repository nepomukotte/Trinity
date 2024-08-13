/* =====================================================================================================================
 * Logging.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class utl::Logger and macros for logging helpers.
 * Created on: Apr 4, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Frederic Chateau (frederic.chateau@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the the Mordicus Real-Time Framework general utilities.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
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
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#include "Logging.h"
#include <string>
#include <cstdlib>

namespace utl
{
// =============================================================================
// Methods of class BackendLogger
// -----------------------------------------------------------------------------

/**
 * The logging back-end.
 * We must use a raw pointer because we are performing lazy initialization, and
 * this is the only way to be sure it is initialized at 0 before any object
 * construction.
 */
AbstractLoggingBackend* BackendLogger::backend_ = 0;

AbstractLoggingBackend* PlainBackendLogger::backend()
{
	if(backend_ == 0)
	{
		setBackend(LoggingBackendPtr(new StdoutLoggingBackend));
	}
	return backend_;
}

void PlainBackendLogger::setBackend(LoggingBackendPtr b)
{
	// only registering cleanup function if backend is null
	// ie: first use, or already destroyed
	// (otherwise we could register it several times)
	if(backend_ == 0)
	{
		// using atexit to implement "phoenix singleton" pattern
		// (see Modern C++ Design by Andrei Alexandrescu)
		std::atexit(destroyBackend);
	}
	else
	{
		destroyBackend();
	}
	backend_ = b.release();
}

void PlainBackendLogger::destroyBackend()
{
	delete backend_;
	backend_ = 0;
}

// =============================================================================
// Methods of class AbstractLoggingBackend
// -----------------------------------------------------------------------------

AbstractLoggingBackend::~AbstractLoggingBackend()
{
}

/**
 * @fn std::ostream& AbstractLoggingBackend::beginLog(SeverityLevel level, const char* source, unsigned line)
 * This method is called whenever a logging expression begins. It shall return a
 * stream where the expression will be appended.
 */

/**
 * @fn void AbstractLoggingBackend::endLog()
 * This method is called whenever a logging expression ends. It shall do
 * something with the stream it previously returned so that the log is
 * displayed/stored/sent or whatever...
 */

// =============================================================================
// Methods of class StdoutLoggingBackend
// -----------------------------------------------------------------------------

std::ostream& StdoutLoggingBackend::beginLog(SeverityLevel /*level*/, const char* /*source*/, unsigned /*line*/)
{
	return std::cout;
}

void StdoutLoggingBackend::endLog()
{
	std::cout << std::endl;
}

// =============================================================================
// Free-standing functions
// -----------------------------------------------------------------------------

std::string xtract_base(const std::string& filePath)
{
	size_t head = filePath.find_last_of('/');
	head = (head == std::string::npos) ? 0u : head+1;
	size_t tail = filePath.find_last_of('.');
	return filePath.substr(head, tail-head);
}

}
