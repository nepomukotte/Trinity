#ifndef CCFG_IO_DB_EXCEPTION_H_INCLUDED
#define CCFG_IO_DB_EXCEPTION_H_INCLUDED

#include "CCfg/Exception.h"

namespace CCfg
{
namespace Io
{
namespace Db
{

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is the base class for exception declared in the database module.
 * @author fchateau
 */
class Exception : public CCfg::Exception
{
public:
	Exception();
	Exception(const std::string& message);
	Exception(const Exception& value);
	Exception& operator =(const Exception& value);
	virtual ~Exception() throw();

protected:
	virtual std::string formatMessage() const = 0;
};

}
}
}

#endif // CCFG_IO_DB_EXCEPTION_H_INCLUDED
