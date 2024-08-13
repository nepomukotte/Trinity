// =============================================================================
// UnsupportedOperation.h
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

#ifndef CCfg_UnsupportedOperation_h_INCLUDED
#define CCfg_UnsupportedOperation_h_INCLUDED

#include "Exception.h"
#include <string>
#include <typeinfo>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This exception class reports that a virtual operation has been called on an
 * object that doesn't support this operation.
 * It makes only sense for class hierarchies where the base class interface is
 * the union of its derived class interfaces instead of the intersection (which is more
 * common).
 *
 * @author fchateau
 */
class UnsupportedOperation : public Exception
{
	std::string name;               ///< the name of the object that doesn't support the operation.
	const std::type_info& clazz;	///< the class of the object.
	std::string operationName;		///< the operation name.
	std::string comment;			///< a comment about the error.

public:
	// construction/destruction
	UnsupportedOperation(const std::string& name,
		const std::type_info& clazz,
		const std::string& operationName,
		const std::string& comment = "");
	UnsupportedOperation(const UnsupportedOperation& value);
	virtual ~UnsupportedOperation() throw();

	// getters
	const std::string& getName() const;
	const std::string& getOperationName() const;
	const std::string& getComment() const;
	const std::type_info& getClass() const;

private:
	std::string formatMessage() const;
	// assignment is not supported because type_info does not support assignment
	UnsupportedOperation& operator =(const UnsupportedOperation& value);
};


///////////////////////////////////////////////////////////////////////////////
// inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Returns \copydoc #name
 */
inline const std::string& UnsupportedOperation::getName() const
{
	return name;
}

/**
 * Getter of the property "operationName".
 * @return Returns the name of the unsupported operation.
 */
inline const std::string& UnsupportedOperation::getOperationName() const
{
	return operationName;
}

/**
 * Getter of the property "comment"
 * @return Returns the comment about the error.
 */
inline const std::string& UnsupportedOperation::getComment() const
{
	return comment;
}

/**
 * Getter of property "clazz"
 * @return Returns the type information about the class reporting the error.
 */
inline const std::type_info& UnsupportedOperation::getClass() const
{
	return clazz;
}

}


#endif // CCfg_UnsupportedOperation_h_INCLUDED
