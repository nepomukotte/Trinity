// =============================================================================
// Comment.h
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

#ifndef CCfg_Comment_h_INCLUDED
#define CCfg_Comment_h_INCLUDED

#include <string>
#include <ostream>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class describes a comment in a configuration document.
 *
 * @author fchateau
 */
class Comment
{
	std::string text;

public:
	// construction/destruction
	Comment();
	Comment(const std::string& commentText);
	Comment(const Comment& value);
	Comment& operator =(const Comment& value);
	~Comment();

	// comparison operators
	bool operator==(const Comment& right) const;

	// getters/setters
	const std::string& getText() const;
	void setText(const std::string& value);
};

std::ostream& operator<<(std::ostream& stream, const Comment& comment);


///////////////////////////////////////////////////////////////////////////////
// inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Getter of the property "text".
 * @return Returns the text of the comment.
 */
inline const std::string& Comment::getText() const
{
	return text;
}

/**
 * Setter of the property "text".
 * @param value the new text for the comment.
 */
inline void Comment::setText(const std::string& value)
{
	text = value;
}

inline std::ostream& operator<<(std::ostream& stream, const Comment& comment)
{
	return stream << comment.getText();
}

}

#endif // CCfg_Comment_h_INCLUDED
