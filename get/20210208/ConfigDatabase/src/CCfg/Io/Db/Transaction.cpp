#include "Transaction.h"
#include <soci.h>
#include <exception>
#include <iostream>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::Transaction class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor: begins a transaction on the specified session.
 */
Transaction::Transaction(SOCI::Session& session)
	: session(session), finished(false)
{
	session.begin();
}

/**
 * Destructor: rolls back unfinished transactions.
 */
Transaction::~Transaction()
{
	if(not finished)
	{
		rollback();
	}
}

/**
 * Commits the transaction.
 */
void Transaction::commit()
{
	session.commit();
	finished = true;
}

/**
 * Rolls back the transaction.
 * Never throws exceptions.
 */
void Transaction::rollback()
{
	try
	{
		session.rollback();
		finished = true;
	}
	catch(std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
	catch(...)
	{
		std::cerr << "unknown exception during rollback !" << std::endl;
	}
}

}
}
}
