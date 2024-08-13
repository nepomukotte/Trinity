#ifndef UTILS_ERROR_ABSTRACTFORMATTER_H_INCLUDED
#define UTILS_ERROR_ABSTRACTFORMATTER_H_INCLUDED

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
 * This base class extracts informations from exception classes and format
 * messages of different semantic nature about the error.
 *
 * @author fchateau
 */
class AbstractFormatter
{
public:
	AbstractFormatter() {}
	virtual ~AbstractFormatter() {}

	virtual std::string formatTitle(std::exception& e);
	virtual std::string formatBrief(std::exception& e);
	virtual std::string formatDetails(std::exception& e);
	virtual std::string formatSuggestion(std::exception& e);

private:
	AbstractFormatter(const AbstractFormatter& rvalue);
	AbstractFormatter& operator=(const AbstractFormatter& rvalue);
};

}
}

#endif // UTILS_ERROR_ABSTRACTFORMATTER_H_INCLUDED
