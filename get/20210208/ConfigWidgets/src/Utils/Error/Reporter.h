#ifndef UTILS_ERROR_REPORTER_H_INCLUDED
#define UTILS_ERROR_REPORTER_H_INCLUDED

#include "AbstractFormatter.h"
#include "AbstractPrinter.h"
#include <memory>

namespace Utils
{
namespace Error
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class centralizes error reporting.
 * Each time an error is detected, whatever its severity, whether it should be
 * logged or pop up a dialog box to the user, it should be reported using this class.
 * The formatting and printing policies can be customized to implement the error
 * reporting. They can be implemented to have a different behaviour according
 * to the exception type using RTTI.
 *
 * This customization is designed to be realized at application scope, in order
 * to implement the application error message handling policy.
 *
 * The informations about the error are transmitted using classes deriving from
 * std::exception. Netherless, it does not force the client to actually throw
 * exceptions. It is possible to instanciate the error class without throwing it,
 * just to transmit the informations, and report it immediately.
 *
 * @author fchateau
 */
class Reporter
{
	static std::auto_ptr<AbstractFormatter> formatter;
	static std::auto_ptr<AbstractPrinter> printer;
public:
	static void report(std::exception& e);
	static void setFormatter(std::auto_ptr<AbstractFormatter> newFormatter);
	static void setPrinter(std::auto_ptr<AbstractPrinter> newPrinter);

private:
	Reporter() {}
	Reporter(const Reporter& rvalue);
	Reporter& operator=(const Reporter& rvalue);
};

}
}

#endif // UTILS_ERROR_REPORTER_H_INCLUDED
