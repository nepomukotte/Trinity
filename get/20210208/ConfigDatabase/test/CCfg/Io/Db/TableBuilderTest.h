#ifndef CCFG_IO_DB_TABLEBUILDERTEST_H_INCLUDED
#define CCFG_IO_DB_TABLEBUILDERTEST_H_INCLUDED

#include "CCfg/Io/Db/TableBuilder.h"
#include <string>

namespace CCfg
{
namespace Io
{
namespace Db
{
std::string createQuery1(TableBuilder& builder);
std::string createQuery2(TableBuilder& builder);
std::string createQuery3(TableBuilder& builder);
std::string createQuery4(TableBuilder& builder);
}
}
}

#endif // CCFG_IO_DB_TABLEBUILDERTEST_H_INCLUDED
