#include "AbstractPrinter.h"
#include <iostream>

namespace Utils
{
namespace Error
{
////////////////////////////////////////////////////////////////////////////////
// Utils::Error::AbstractPrinter class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * This methods implements the printing of the error message, using the formatter
 * to produce error description messages.
 * The default implementation prints the brief description and the detailed
 * description on stderr.
 */
void AbstractPrinter::print(std::exception& e, AbstractFormatter& formatter)
{
	std::string brief = formatter.formatBrief(e);
	std::string details = formatter.formatDetails(e);
	if(not brief.empty())
	{
		std::cerr << brief << std::endl;
	}
	if(not details.empty())
	{
		std::cerr << details << std::endl;
	}
}

}
}
