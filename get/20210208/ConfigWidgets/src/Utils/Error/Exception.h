#ifndef UTILS_ERROR_EXCEPTION_H_INCLUDED
#define UTILS_ERROR_EXCEPTION_H_INCLUDED

#include <string>
#include <exception>

namespace Utils
{
namespace Error
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is a base exception class providing methods to describe the error more
 * accurately than a classic exception class.
 *
 * @author fchateau
 */
class Exception : public std::exception
{
	mutable std::string title;
	mutable std::string brief;
	mutable std::string details;
	mutable std::string suggestion;
	mutable std::string message;
public:
	Exception();
	Exception(const Exception& rvalue);
	virtual ~Exception() throw();
	Exception& operator=(const Exception& rvalue);

	const char* what() const throw();

	const std::string& getTitle() const;
	const std::string& getBrief() const;
	const std::string& getDetails() const;
	const std::string& getSuggestion() const;
	const std::string& getMessage() const;

protected:
	virtual std::string formatTitle() const;
	virtual std::string formatBrief() const = 0;
	virtual std::string formatDetails() const;
	virtual std::string formatSuggestion() const;
	virtual std::string formatMessage() const;
};

}
}

#endif // UTILS_ERROR_EXCEPTION_H_INCLUDED
