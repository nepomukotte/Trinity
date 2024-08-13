/* =====================================================================================================================
 * Logging.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class utl::Logger and macros for logging helpers.
 * Created on: Jan 11, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef utl_Logging_h_INCLUDED
#define utl_Logging_h_INCLUDED

#ifdef UTL_LOG_THREADSAFE
#	include <utl/Singleton.hpp>
#	include <utl/Singleton.ipp>
#	ifdef VXWORKS
#		include <vxWorks.h>
#		include <semLib.h>
#	else
#		include <boost/thread.hpp>
#	endif /* VXWORKS */
#endif /* UTL_LOG_THREADSAFE */


#include <string>
#include <sstream>
#include <iostream>
#include <memory>

namespace utl
{
enum SeverityLevel
{
	debug = 0,
	info,
	warn,
	error,
	fatal
};

class AbstractLoggingBackend;

typedef std::auto_ptr<AbstractLoggingBackend> LoggingBackendPtr;

/**
 * This class allows to declare log expressions that will be processed by
 * a parameterizable global backend.
 * It is supposed to be instantiated at the beginning of each expression and
 * destroyed at the end. It is not thread-safe
 */
class PlainBackendLogger
{
public:
	PlainBackendLogger(SeverityLevel level, const char* source, unsigned line);
	~PlainBackendLogger();

	template<class T>
	PlainBackendLogger& operator<<(const T& val);

	static AbstractLoggingBackend* backend();
	static void setBackend(LoggingBackendPtr b);

protected:
	static void destroyBackend();

	static AbstractLoggingBackend* backend_;
	std::ostream& stream;
};

#ifdef UTL_LOG_THREADSAFE

#ifndef VXWORKS
/**
 * This structure embeds a boost mutex so as to make it possible to instantiate the mutex in a utl::Singleton template.
 */
struct LoggingMutex
{
	boost::recursive_mutex instance;
};

typedef boost::unique_lock<boost::recursive_mutex> MutexLock;

#define LOGGING_MUTEX ::utl::Singleton<LoggingMutex>::instance().instance

#else // i.e. VxWorks context
/**
 * This structure embeds a boost mutex so as to make it possible to instantiate the mutex in a utl::Singleton template.
 */
class LoggingMutex
{
	friend class MutexLock;

public:
	LoggingMutex() : vxInstance(semMCreate(SEM_Q_FIFO)) {} // All mutexes are recursive in VxWorks

private:
	SEM_ID vxInstance;
};

class MutexLock
{
public:
	MutexLock(LoggingMutex& mtx) : mtx(mtx)
	{
		semTake(mtx.vxInstance, WAIT_FOREVER);
	}

	~MutexLock()
	{
		semGive(mtx.vxInstance);
	}

private:

	LoggingMutex& mtx;
};

#define LOGGING_MUTEX ::utl::Singleton<LoggingMutex>::instance()

#endif /* VXWORKS */


/**
 * This class is a thread-safe version of BackendLogger.
 */
class ThreadSafeBackendLogger
{
public:
	ThreadSafeBackendLogger(SeverityLevel level, const char* source, unsigned line)
	: nowLogging(LOGGING_MUTEX), backendLogger(level, source, line) {}

	template<class T>
	PlainBackendLogger& operator<<(const T& val)
	{
		return backendLogger << val;
	}

	static AbstractLoggingBackend* backend()
	{
		MutexLock l(LOGGING_MUTEX);
		return PlainBackendLogger::backend();
	}

	static void setBackend(LoggingBackendPtr b)
	{
		MutexLock l(LOGGING_MUTEX);
		PlainBackendLogger::setBackend(b);
	}

private:
	MutexLock nowLogging;
	PlainBackendLogger backendLogger;
};

/**
 * In case of thread-safe logging, the back-end logger is defined as the more complex tread-safe class.
 */
typedef ThreadSafeBackendLogger BackendLogger;
#undef LOGGING_MUTEX

#else /* UTL_LOG_THREADSAFE */
/**
 * In case of non thread-safe logging, the back-end logger is defined as the plain non tread-safe class.
 */
typedef PlainBackendLogger BackendLogger;

#endif /* UTL_LOG_THREADSAFE */

/**
 * This class implements a very simple log expression that just outputs the
 * expression on stdout and appends it an endline.
 */
class StdoutLogger
{
public:
	~StdoutLogger() { std::cout << std::endl; }
	template<class T>
	StdoutLogger& operator<<(const T& val)
	{
		std::cout << val;
		return *this;
	}
};

/**
 * This class defines the interface for logging backends.
 * Classes deriving from it will have to implement how the log lines shall be
 * displayed/stored.
 */
class AbstractLoggingBackend
{
public:
	virtual ~AbstractLoggingBackend() = 0;
	virtual std::ostream& beginLog(SeverityLevel level, const char* source, unsigned line) = 0;
	virtual void endLog() = 0;
};

/**
 * This is a logging backend based on stdout. It ensures logs are written on stdout.
 * It is used as the default backend implementation for BackendLogger.
 */
class StdoutLoggingBackend : public AbstractLoggingBackend
{
public:
	std::ostream& beginLog(SeverityLevel level, const char* source, unsigned line);
	void endLog();
};

std::string xtract_base(const std::string& filePath);

/**
 * This class bufferizes a log expression and only sends it when the object
 * goes out of scope.
 */
class ScopeLog
{
public:
	ScopeLog(SeverityLevel level, const char* source, unsigned line)
		: level_(level), source_(source), line_(line)
	{
	}
	~ScopeLog()
	{
		PlainBackendLogger logger(level_, source_, line_);
		logger << stream_.str();
	}
	std::ostream& stream() { return stream_; }
private:
	SeverityLevel level_;
	const char* source_;
	unsigned line_;
	std::ostringstream stream_;
};

// =============================================================================
// Inline methods of class BackendLogger
// -----------------------------------------------------------------------------
#ifndef __CINT__ // Let the ROOT dictionary generator rootcint  skip this
inline PlainBackendLogger::PlainBackendLogger(SeverityLevel level, const char* source, unsigned line)
	: stream(backend()->beginLog(level, source, line))
{
}

inline PlainBackendLogger::~PlainBackendLogger()
{
	backend_->endLog();
}

template<class T>
inline PlainBackendLogger& PlainBackendLogger::operator<<(const T& val)
{
	stream << val;
	return *this;
}
#endif // __CINT__
}

// =============================================================================
// General Macros
// -----------------------------------------------------------------------------

#ifndef UTL_LOG_LEVEL
# ifdef UTL_LOGGING_ENABLED
#  define UTL_LOG_LEVEL 0
# else
#  define UTL_LOG_LEVEL 10
# endif
#else
# ifndef UTL_LOGGING_ENABLED
#  define UTL_LOGGING_ENABLED
# endif
#endif


// =============================================================================
// LOG Macro
// -----------------------------------------------------------------------------

#ifdef UTL_LOGGING_ENABLED
# define LOG(expr) ::utl::StdoutLogger() << expr;
#else
# define LOG(expr)
#endif

// =============================================================================
// LOG_ Macros:
// Simple macros to format and send immediately a log line according to its
// severity level.
// It relies on branch optimization to evict disabled log messages. Some compilers
// may not optimize it away when the log level condition is not met.
// -----------------------------------------------------------------------------

#define LOG_DEBUG() if(UTL_LOG_LEVEL > 0) ; else ::utl::BackendLogger(::utl::debug, __FILE__, __LINE__)
#define LOG_INFO()  if(UTL_LOG_LEVEL > 1) ; else ::utl::BackendLogger(::utl::info,  __FILE__, __LINE__)
#define LOG_WARN()  if(UTL_LOG_LEVEL > 2) ; else ::utl::BackendLogger(::utl::warn,  __FILE__, __LINE__)
#define LOG_ERROR() if(UTL_LOG_LEVEL > 3) ; else ::utl::BackendLogger(::utl::error, __FILE__, __LINE__)
#define LOG_FATAL() if(UTL_LOG_LEVEL > 4) ; else ::utl::BackendLogger(::utl::fatal, __FILE__, __LINE__)

// =============================================================================
// MSG_ Macros
// These macros format and send a log line whatever the compile-time log level is.
// -----------------------------------------------------------------------------

#define MSG() ::utl::BackendLogger(::utl::info,  __FILE__, __LINE__)
#define MSG_IF(cond)  if(!(cond)) ; else ::utl::BackendLogger(::utl::info,  __FILE__, __LINE__)

// =============================================================================
// MLOG_ Macros
// Same as LOG_ macros, except the need to put the whole expression in
// parentheses, and the fact it relies on macro evaluation rather than branch
// optimization. So it is guaranteed to evaluate to nothing on any compiler if
// the compile-time log level is higher.
// It also enables referencing symbols that would be deactivated by the
// preprocessor on higher log levels (eg: a member variable surrounded by
// a #if UTL_LOG_LEVEL < 1), resulting in a compilation error with LOG_ macros.
// -----------------------------------------------------------------------------

#if UTL_LOG_LEVEL < 1
# define MLOG_DEBUG(expr) ::utl::BackendLogger(::utl::debug, __FILE__, __LINE__) << expr;
#endif
#if UTL_LOG_LEVEL < 2
# define MLOG_INFO(expr)  ::utl::BackendLogger(::utl::info,  __FILE__, __LINE__) << expr;
#endif
#if UTL_LOG_LEVEL < 3
# define MLOG_WARN(expr)  ::utl::BackendLogger(::utl::warn,  __FILE__, __LINE__) << expr;
#endif
#if UTL_LOG_LEVEL < 4
# define MLOG_ERROR(expr) ::utl::BackendLogger(::utl::error, __FILE__, __LINE__) << expr;
#endif
#if UTL_LOG_LEVEL < 5
# define MLOG_FATAL(expr) ::utl::BackendLogger(::utl::fatal, __FILE__, __LINE__) << expr;
#endif

#ifndef MLOG_DEBUG
# define MLOG_DEBUG(expr)
#endif
#ifndef MLOG_INFO
# define MLOG_INFO(expr)
#endif
#ifndef MLOG_WARN
# define MLOG_WARN(expr)
#endif
#ifndef MLOG_ERROR
# define MLOG_ERROR(expr)
#endif
#ifndef MLOG_FATAL
# define MLOG_FATAL(expr)
#endif

// =============================================================================
// SLOG_* Macros
// Scope log macro: the SLOG_*_BEGIN macros initialize a log line which will be
// sent at the end of the scope even in case of exception.
// It is possible to concatenate the log line using SLOG_ macros until the end
// of the scope.
// -----------------------------------------------------------------------------

#if UTL_LOG_LEVEL < 1
# define SLOG_DEBUG_BEGIN(var) ::utl::ScopeLog var(::utl::debug, __FILE__, __LINE__); var.stream()
#endif
#if UTL_LOG_LEVEL < 2
# define SLOG_INFO_BEGIN(var)  ::utl::ScopeLog var(::utl::info,  __FILE__, __LINE__); var.stream()
#endif
#if UTL_LOG_LEVEL < 3
# define SLOG_WARN_BEGIN(var)  ::utl::ScopeLog var(::utl::warn,  __FILE__, __LINE__); var.stream()
#endif
#if UTL_LOG_LEVEL < 4
# define SLOG_ERROR_BEGIN(var) ::utl::ScopeLog var(::utl::error, __FILE__, __LINE__); var.stream()
#endif
#if UTL_LOG_LEVEL < 5
# define SLOG_FATAL_BEGIN(var) ::utl::ScopeLog var(::utl::fatal, __FILE__, __LINE__); var.stream()
#endif

#ifndef SLOG_DEBUG_BEGIN
# define SLOG_DEBUG_BEGIN(var) LOG_DEBUG()
#endif
#ifndef SLOG_INFO_BEGIN
# define SLOG_INFO_BEGIN(var)  LOG_INFO()
#endif
#ifndef SLOG_WARN_BEGIN
# define SLOG_WARN_BEGIN(var)  LOG_WARN()
#endif
#ifndef SLOG_ERROR_BEGIN
# define SLOG_ERROR_BEGIN(var) LOG_ERROR()
#endif
#ifndef SLOG_FATAL_BEGIN
# define SLOG_FATAL_BEGIN(var) LOG_FATAL()
#endif

#if UTL_LOG_LEVEL < 1
# define SLOG_DEBUG(var) var.stream()
#endif
#if UTL_LOG_LEVEL < 2
# define SLOG_INFO(var)  var.stream()
#endif
#if UTL_LOG_LEVEL < 3
# define SLOG_WARN(var)  var.stream()
#endif
#if UTL_LOG_LEVEL < 4
# define SLOG_ERROR(var) var.stream()
#endif
#if UTL_LOG_LEVEL < 5
# define SLOG_FATAL(var) var.stream()
#endif

#ifndef SLOG_DEBUG
# define SLOG_DEBUG(var) LOG_DEBUG()
#endif
#ifndef SLOG_INFO
# define SLOG_INFO(var)  LOG_INFO()
#endif
#ifndef SLOG_WARN
# define SLOG_WARN(var)  LOG_WARN()
#endif
#ifndef SLOG_ERROR
# define SLOG_ERROR(var) LOG_ERROR()
#endif
#ifndef SLOG_FATAL
# define SLOG_FATAL(var) LOG_FATAL()
#endif

// =============================================================================
// PROMPT_ Macros
// -----------------------------------------------------------------------------

#define PROMPT_MSG(expr) LOG("--- " << ::utl::xtract_base(__FILE__) << ": " << expr)
#define PROMPT_ERR(expr) LOG("*** " << ::utl::xtract_base(__FILE__) << ": " << expr)

#endif // utl_Logging_h_INCLUDED
