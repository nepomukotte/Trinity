// =============================================================================
// Exception.h
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef CCfg_Exception_h_INCLUDED
#define CCfg_Exception_h_INCLUDED

#include <string>
#include <exception>

namespace CCfg
{

/**
 * This class is the base class of every exception in the module.
 * It must be derived in order to define specific error informations
 * for each exception type. The derived class must implement the formatMessage
 * method because it will be called used to initialize the message attribute.
 *
 * @author fchateau
 */
class Exception : public std::exception
{
	mutable std::string message;	///< the error message. "mutable" is for lazy initialization in getter.

public:
	Exception();
	explicit Exception(const std::string& message);
	Exception(const Exception& value);
	Exception& operator =(const Exception& value);
	virtual ~Exception() throw();

	const char* what() const throw();

	const std::string& getMessage() const;

protected:
	/**
	 * Format an error message for the actual exception based on its attributes.
	 * @return Returns the generated error message.
	 */
	virtual std::string formatMessage() const = 0;
};






///////////////////////////////////////////////////////////////////////////////
// inline methods body
///////////////////////////////////////////////////////////////////////////////

}


#endif // CCfg_Exception_h_INCLUDED
