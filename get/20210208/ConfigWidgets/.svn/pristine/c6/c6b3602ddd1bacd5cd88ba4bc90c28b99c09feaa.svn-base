#ifndef UTILS_ERROR_ABSTRACTPRINTER_H_INCLUDED
#define UTILS_ERROR_ABSTRACTPRINTER_H_INCLUDED

#include "AbstractFormatter.h"

namespace Utils
{
namespace Error
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This base class can define how/where/which error messages produced by the
 * formatter are printed/displayed.
 *
 * @author fchateau
 */
class AbstractPrinter
{
public:
	AbstractPrinter() {}
	virtual ~AbstractPrinter() {}

	virtual void print(std::exception& e, AbstractFormatter& formatter);

private:
	AbstractPrinter(const AbstractPrinter& rvalue);
	AbstractPrinter& operator=(const AbstractPrinter& rvalue);
};

}
}

#endif // UTILS_ERROR_ABSTRACTPRINTER_H_INCLUDED
