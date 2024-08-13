// =============================================================================
// ServantDescriptor.h
// -----------------------------------------------------------------------------
// class dhsm::ServantDescriptor
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_ServantDescriptor_h_INCLUDED
#define dhsm_ServantDescriptor_h_INCLUDED

#include <string>
#include <memory>

namespace dhsm
{
/**
 * This class describes a servant state-machine engine managed by a MasterEngine.
 */
class ServantDescriptor
{
public:
	class Identity;

	enum Status
	{
		Working = 0,  ///< the servant is up and running. The state machine is in sync with the master.
		Failing,      ///< the servant engine has failed to execute a previous command, and is now out of sync.
		Critical,     ///< the servant engine has raised a critical error
		NotResponding,///< the servant engine has not responded to the previous command, it may be down.
		Unknown       ///< the servant is in an unknown state
	};

	ServantDescriptor();
	virtual ~ServantDescriptor();

	Status status() const { return _status; }
	void setStatus(Status newStatus) { _status = newStatus; }

	const std::string& errorMessage() const { return _errorMessage; }
	void setErrorMessage(const std::string& msg) { _errorMessage = msg; }

	bool isOk() const { return _status == Working; }
	bool isBad() const { return _status == Failing; }
	bool isLost() const { return _status == NotResponding or _status == Critical or _status == Unknown; }

	virtual const Identity& identity() const = 0;

private:
	// non-copyable
	ServantDescriptor(const ServantDescriptor& rvalue);
	ServantDescriptor& operator=(const ServantDescriptor& rvalue);

	Status _status;				///< the servant engine status
	std::string _errorMessage;	///< the error message if the servant is lost
};

class ServantDescriptor::Identity
{
public:
	virtual ~Identity();
	friend std::ostream& operator<<(std::ostream& stream, const Identity& servantId);

protected:
	virtual void write(std::ostream& stream) const = 0;
};

typedef std::auto_ptr<ServantDescriptor> ServantDescriptorPtr;

}

#endif // dhsm_ServantDescriptor_h_INCLUDED
