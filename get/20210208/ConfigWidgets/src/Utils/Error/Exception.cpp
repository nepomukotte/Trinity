#include "Exception.h"

namespace Utils
{
namespace Error
{
////////////////////////////////////////////////////////////////////////////////
// Utils::Error::Exception class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
Exception::Exception()
{
}

/**
 * Copy constructor.
 */
Exception::Exception(const Exception& rvalue)
	: std::exception(rvalue), title(rvalue.title), brief(rvalue.brief), details(rvalue.details),
	suggestion(rvalue.suggestion), message(rvalue.message)
{
}

/**
 * Destructor.
 */
Exception::~Exception() throw()
{
}

/**
 * Assignment operator.
 */
Exception& Exception::operator=(const Exception& rvalue)
{
	title = rvalue.title;
	brief = rvalue.brief;
	details = rvalue.details;
	suggestion = rvalue.suggestion;
	message = rvalue.message;
	return *this;
}

const char* Exception::what() const throw()
{
	return getMessage().c_str();
}

/**
 * Returns the title of the error.
 */
const std::string& Exception::getTitle() const
{
	if(title.empty())
	{
		title = formatTitle();
	}
	return title;
}

/**
 * Returns a brief (one line) description of the error.
 */
const std::string& Exception::getBrief() const
{
	if(brief.empty())
	{
		brief = formatBrief();
	}
	return brief;
}

/**
 * Returns a detailed description of the error.
 */
const std::string& Exception::getDetails() const
{
	if(details.empty())
	{
		details = formatDetails();
	}
	return details;
}

/**
 * Returns a suggestion message to solve the error.
 * The suggestion is destined to the end-user.
 */
const std::string& Exception::getSuggestion() const
{
	if(suggestion.empty())
	{
		suggestion = formatSuggestion();
	}
	return suggestion;
}

/**
 * Returns the standard error message used by what().
 * By default, it is the brief error message.
 */
const std::string& Exception::getMessage() const
{
	if(message.empty())
	{
		message = formatMessage();
	}
	return message;
}

/**
 * Generates a title for this error.
 * Default implementation returns an empty string.
 */
std::string Exception::formatTitle() const
{
	return std::string();
}

/**
 * @fn std::string Exception::formatBrief() const
 * Generates a brief description for this error.
 */

/**
 * Generates a detailed description for this error.
 * Default implementation returns an empty string.
 */
std::string Exception::formatDetails() const
{
	return std::string();
}

/**
 * Generates a suggestion message to resolve this error.
 * Default implementation returns an empty string.
 */
std::string Exception::formatSuggestion() const
{
	return std::string();
}

/**
 * Generates an regular error message to be returned by what().
 * Default implementation uses the brief error message.
 */
std::string Exception::formatMessage() const
{
	return formatBrief();
}

}
}
