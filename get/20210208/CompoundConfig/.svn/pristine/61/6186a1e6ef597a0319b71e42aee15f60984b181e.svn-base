// =============================================================================
// Comment.cpp
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

#include "Comment.h"

namespace CCfg
{

/**
 * Default constructor.
 */
Comment::Comment()
{
}

/**
 * Implicit constructor: constructs a comment with the specified text.
 * @param commentText the text of the comment.
 */
Comment::Comment(const std::string& commentText)
	: text(commentText)
{
}

/**
 * Copy constructor.
 * @param value the comment to copy.
 */
Comment::Comment(const Comment& value)
	: text(value.text)
{
}

/**
 * Destructor.
 */
Comment::~Comment()
{
}

/**
 * Assignment operator.
 * @param value the comment to copy.
 */
Comment& Comment::operator =(const Comment& value)
{
	setText(value.getText());
	return *this;
}

/**
 * Equality operator.
 * @return Returns true if comments are equals, false otherwise.
 */
bool Comment::operator==(const Comment& right) const
{
	return text == right.text;
}


}

