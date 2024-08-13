#ifndef CCFG_IO_DB_TRANSACTION_H_INCLUDED
#define CCFG_IO_DB_TRANSACTION_H_INCLUDED

// forward declarations
namespace SOCI { class Session; }

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
 * This is a wrapper to have exception safe RAII compliant transactions.
 *
 * @author fchateau
 */
class Transaction
{
	SOCI::Session& session;
	bool finished;

public:
	Transaction(SOCI::Session& session);
	virtual ~Transaction();

	void commit();
	void rollback();

private:
	Transaction(const Transaction& rvalue);
	Transaction& operator=(const Transaction& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_TRANSACTION_H_INCLUDED
