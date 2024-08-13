#include "Reporter.h"

namespace Utils
{
namespace Error
{
////////////////////////////////////////////////////////////////////////////////
// Utils::Error::Reporter class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * This object is responsible for extracting a meaningfull and formatted message
 * from exceptions.
 */
std::auto_ptr<AbstractFormatter> Reporter::formatter;
/**
 * This object is responsible for displaying the formatted message to the user,
 * or printing it somewhere.
 */
std::auto_ptr<AbstractPrinter> Reporter::printer;

/**
 * Reports the specified exception, in an appropriate way according to the
 * nature of the exception.
 * The exact behaviour of this method depends on the class of the formatter and
 * printer objects.
 * The default for exceptions deriving from Utils::Error::Exception is to print
 * the brief and detailed messages on stderr.
 * For other exceptions it prints what() on stderr.
 */
void Reporter::report(std::exception& e)
{
	if(formatter.get() == 0)
	{
		formatter.reset(new AbstractFormatter);
	}
	if(printer.get() == 0)
	{
		printer.reset(new AbstractPrinter);
	}
	printer->print(e, *formatter);
}

/**
 * Changes the formatting policy for reported errors.
 */
void Reporter::setFormatter(std::auto_ptr<AbstractFormatter> newFormatter)
{
	formatter = newFormatter;
}

/**
 * Changes the printing policy for reported errrors.
 */
void Reporter::setPrinter(std::auto_ptr<AbstractPrinter> newPrinter)
{
	printer = newPrinter;
}

}
}
