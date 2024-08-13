// =============================================================================
// ValidationVisitor.cpp
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

#include "ValidationVisitor.h"
#include "Object.h"
#include "Scalar.h"
#include "Range.h"
#include "View/ClassIterator.h"
#include <sstream>

namespace CCfg
{

////////////////////////////////////////////////////////////////////////////////
// CCfg::AbstractValidationVisitor class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
AbstractValidationVisitor::AbstractValidationVisitor(int checkFlags)
	: checkFlags(checkFlags), checkedClasses()
{
}

/**
 * Destructor
 */
AbstractValidationVisitor::~AbstractValidationVisitor()
{
}

/*
 *
 */
bool AbstractValidationVisitor::onObjectEnter(const Object& objectCfg)
{
	if(checkFlags & CheckNoDoublets)
	{
		checkDoublet(objectCfg);
	}
	return true;
}

/*
 *
 */
bool AbstractValidationVisitor::onObjectExit(const Object& /*objectCfg*/)
{
	return true;
}

/*
 *
 */
bool AbstractValidationVisitor::onScalar(const Scalar& scalar)
{
	if(checkFlags & CheckNoDoublets)
	{
		checkDoublet(scalar);
	}
	if(checkFlags & CheckRanges)
	{
		checkRange(scalar);
	}
	if(checkFlags & CheckTypes)
	{
		checkType(scalar);
	}
	return true;
}

/**
 * Check if there are doublets of the specified attribute, and report the problem.
 */
void AbstractValidationVisitor::checkDoublet(const Attribute& attribute)
{
	Object* parent = attribute.getParent();
	if(parent != 0)
	{
		for(Attribute::const_iterator sibling = parent->begin(); sibling != parent->end(); ++sibling)
		{
			if(sibling != &attribute and
				sibling->getName() == attribute.getName() and
				sibling->getIndex() == attribute.getIndex())
			{
				onDoublet(attribute); // notify that the checked attribute has doublet(s)
				break; // no need to search more
			}
		}
	}
}

/**
 * Check if the specified scalar is out of range, and report the problem.
 */
void AbstractValidationVisitor::checkRange(const Scalar& scalar)
{
	const ClassInfoPtr& classInfo = scalar.getClassInfo();
	if(classInfo)
	{
		Range* range = scalar.getRange();
		if(range != 0 and not range->contains(scalar.getValue()))
		{
			onOutOfRange(scalar);
		}
	}
}

/**
 * Check if all the attributes of the equivalent class of the specified attribute have the
 * same type and report the problem.
 */
void AbstractValidationVisitor::checkType(const Attribute& attr)
{
	std::string classPath = attr.getClassPath().toString();
	if(find(checkedClasses.begin(), checkedClasses.end(), classPath) == checkedClasses.end())
	{
		View::ClassIterator classIt(const_cast<Attribute&>(attr));
		while(classIt.hasNext())
		{
			classIt.next();
			if((*classIt.get()).getType() != attr.getType())
			{
				onTypeMismatch(attr, *classIt.get());
			}
		}
		checkedClasses.push_back(classPath);
	}
}

/**
 * @fn void AbstractValidationVisitor::onDoublet(const Attribute& attribute)
 * Callback method called each time the validation process encounters a doublet.
 * It will be called once for each attribute involved.
 * @param attribute an attribute having a doublet. Beware: it is NOT the first
 * doublet found, it is the attribute for which we can find a doublet.
 */

/**
 * @fn void AbstractValidationVisitor::onOutOfRange(const Scalar& scalar)
 * Callback method called each time the validation process encounters a scalar
 * value out of its range.
 */

/**
 * @fn void AbstractValidationVisitor::onTypeMismatch(const Attribute& attr, const Attribute& other)
 * Callback method called each time the validation process encounters an
 * equivalence class where the type is not homogeneous.
 */


////////////////////////////////////////////////////////////////////////////////
// CCfg::ValidationVisitor class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ValidationVisitor::ValidationVisitor(int checkFlags)
	: AbstractValidationVisitor(checkFlags)
{
}

/**
 * Destructor.
 */
ValidationVisitor::~ValidationVisitor()
{
}

/*
 *
 */
void ValidationVisitor::onDoublet(const Attribute& attribute)
{
	throw DoubletError(attribute);
}

/*
 *
 */
void ValidationVisitor::onOutOfRange(const Scalar& scalar)
{
	throw OutOfRange(scalar);
}

/*
 *
 */
void ValidationVisitor::onTypeMismatch(const Attribute& attribute, const Attribute& other)
{
	throw TypeMismatch(attribute, other);
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::ValidationVisitor::DoubletError class methods body
////////////////////////////////////////////////////////////////////////////////

ValidationVisitor::DoubletError::DoubletError(const Attribute& attribute)
	: ValidationError(), doubletPath(attribute.getPath())
{
}

std::string ValidationVisitor::DoubletError::formatMessage() const
{
	std::ostringstream s;
	s << "The attribute \"" << doubletPath << "\" has a doublet." << std::endl;
	return s.str();
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::ValidationVisitor::OutOfRange class methods body
////////////////////////////////////////////////////////////////////////////////

ValidationVisitor::OutOfRange::OutOfRange(const Scalar& scalar)
	: ValidationError(), path(scalar.getPath()), value(scalar.getValue()),
	range(scalar.getRange()->clone())
{
}

ValidationVisitor::OutOfRange::OutOfRange(const OutOfRange& rvalue)
	: ValidationError(rvalue), path(rvalue.path), value(rvalue.value),
	range(rvalue.range->clone())
{
}

std::string ValidationVisitor::OutOfRange::formatMessage() const
{
	std::ostringstream s;
	s << "The value of the scalar \"" << path << "\" is out of range." << std::endl;
	return s.str();
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::ValidationVisitor::TypeMismatch class methods body
////////////////////////////////////////////////////////////////////////////////

ValidationVisitor::TypeMismatch::TypeMismatch(const Attribute& first, const Attribute& second)
	: ValidationError(), firstAttrPath(first.getPath()), secondAttrPath(second.getPath()),
	firstAttrType(first.getType()), secondAttrType(second.getType())
{
}

std::string ValidationVisitor::TypeMismatch::formatMessage() const
{
	std::ostringstream s;
	s << "The types of \"" << firstAttrPath << "\" and \"" << secondAttrPath << "\" do not match.";
	return s.str();
}

}
