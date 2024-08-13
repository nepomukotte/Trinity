// =============================================================================
// ValidationVisitor.h
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

#ifndef CCfg_ValidationVisitor_h_INCLUDED
#define CCfg_ValidationVisitor_h_INCLUDED

#include "Visitor.h"
#include "Path.h"
#include "Range.h"
#include "Value.h"
#include "Exception.h"
#include <vector>

namespace CCfg
{
// forward declarations
class Attribute;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is the base class for visitors that validate a configuration.
 * It needs to be overriden to implement callback methods reporting the
 * validation errors.
 *
 * @author fchateau
 */
class AbstractValidationVisitor : public Visitor
{
public:
	// enumerations
	enum CheckFlags {CheckNoDoublets = 1, CheckRanges = 2, CheckTypes = 4};

	// @{ construction
	AbstractValidationVisitor(int checkFlags = CheckNoDoublets | CheckRanges | CheckTypes);
	virtual ~AbstractValidationVisitor();
	// @}

	// @{ overriden visitor callbacks
	bool onObjectEnter(const Object& objectCfg);
	bool onObjectExit(const Object& objectCfg);
	bool onScalar(const Scalar& scalar);
	// @}

	// @{ callback methods for validation error reporting
	virtual void onDoublet(const Attribute& attribute) = 0;
	virtual void onOutOfRange(const Scalar& scalar) = 0;
	virtual void onTypeMismatch(const Attribute& attr, const Attribute& other) = 0;
	// @}

protected:
	// typedefs
	typedef std::vector<std::string> CheckedClassesList;

	int checkFlags;
	CheckedClassesList checkedClasses;
private:
	// @{ non-copyable
	AbstractValidationVisitor(const AbstractValidationVisitor& rvalue);
	AbstractValidationVisitor& operator=(const AbstractValidationVisitor& rvalue);
	// @}

	// @{ validation methods
	void checkDoublet(const Attribute& attribute);
	void checkRange(const Scalar& scalar);
	void checkType(const Attribute& attribute);
	// @}
};

/**
 * This class provides a straightforward implementation of the AbstractValidationVisitor
 * that just throws an exception when it finds a validation error.
 */
struct ValidationVisitor : public AbstractValidationVisitor
{
	// inner classes
	class ValidationError;
	class DoubletError;
	class OutOfRange;
	class TypeMismatch;

	// @{ construction
	ValidationVisitor(int checkFlags = CheckNoDoublets | CheckRanges | CheckTypes);
	virtual ~ValidationVisitor();
	// @}

	// @{ validation callbacks implementation
	void onDoublet(const Attribute& attribute);
	void onOutOfRange(const Scalar& scalar);
	void onTypeMismatch(const Attribute& attr, const Attribute& other);
	// @}

private:
	ValidationVisitor(const ValidationVisitor& rvalue);
	ValidationVisitor& operator=(const ValidationVisitor& rvalue);
};

/**
 *
 */
class ValidationVisitor::ValidationError : public Exception
{
public:
	ValidationError() : Exception() {}
	virtual ~ValidationError() throw() {}
};

/**
 *
 */
class ValidationVisitor::DoubletError : public ValidationError
{
	Path doubletPath;	///< the path of the attribute that have a doublet.
public:
	DoubletError(const Attribute& attribute);
	~DoubletError() throw() {}

	const Path& getDoubletPath() const { return doubletPath; }

protected:
	std::string formatMessage() const;
};

/**
 *
 */
class ValidationVisitor::OutOfRange : public ValidationError
{
	Path path; ///< the path of the out of range scalar.
	Value value;	///< the out of range value.
	RangePtr range;	///< the range.
public:
	OutOfRange(const Scalar& scalar);
	OutOfRange(const OutOfRange& rvalue);
	~OutOfRange() throw() {}

	const Path& getPath() const { return path; }
	const Value& getValue() const { return value; }
	const Range& getRange() const { return *range; }

protected:
	std::string formatMessage() const;
};

/**
 *
 */
class ValidationVisitor::TypeMismatch : public ValidationError
{
	Path firstAttrPath, secondAttrPath;
	AttributeTypeId firstAttrType, secondAttrType;
public:
	TypeMismatch(const Attribute& first, const Attribute& second);
	~TypeMismatch() throw() {}

	const Path& getFirstPath() const { return firstAttrPath; }
	const Path& getSecondPath() const { return secondAttrPath; }
	AttributeTypeId getFirstType() const { return firstAttrType; }
	AttributeTypeId getSecondType() const { return secondAttrType; }

protected:
	std::string formatMessage() const;
};

}

#endif // CCfg_ValidationVisitor_h_INCLUDED
