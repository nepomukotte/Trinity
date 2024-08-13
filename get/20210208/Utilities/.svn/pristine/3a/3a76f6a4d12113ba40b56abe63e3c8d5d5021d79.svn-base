// =============================================================================
// Exception.h
// -----------------------------------------------------------------------------
// class utl::Exception
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef utl_Exception_INCLUDED
#define utl_Exception_INCLUDED

#include <string>
#include <ostream>
#include <exception>

namespace utl
{
/**
 * This class is the base class of every exception in the module.
 * It must be derived in order to define specific error informations
 * for each exception type. The derived class must implement the writeMessage
 * method because it will be called used to initialize the message attribute.
 *
 * @author fchateau
 */
class Exception : public std::exception
{
	mutable std::string message;	///< the error message. "mutable" is for lazy initialization in getter.

public:
	Exception();
	Exception(const Exception& value);
	Exception& operator =(const Exception& value);
	virtual ~Exception() throw();

	virtual void writeMessage(std::ostream& stream) const = 0;

	const char* what() const throw();
	const std::string& getMessage() const;
};

inline std::ostream& operator<<(std::ostream& stream, const Exception& e)
{
	e.writeMessage(stream);
	return stream;
}

/**
 * Macro for automatic management of retry of a boolean operation.
 * @author Shebli Anvar
 */
#define RETRY(maxRetries, expression, T_Exception, exceptionParam) \
	for(int __cnt=1u; not(expression); __cnt++) {\
		if(__cnt == maxRetries) throw T_Exception(exceptionParam); \
	}

} // namespace utl


#endif // utl_Exception_INCLUDED
