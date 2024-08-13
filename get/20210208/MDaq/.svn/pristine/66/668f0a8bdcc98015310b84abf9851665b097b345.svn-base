/* =============================================================================
 * BaseDynRegister.h
 * -----------------------------------------------------------------------------
 * class mdaq::hw::server::BaseDynRegister
 * class mdaq::hw::server::BaseDynRegister::FieldNotFound
 * Created on: 21 August 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of
 * the Mordicus Real-Time Software Development Framework.
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#ifndef mdaq_hw_server_BaseDynRegister_h_INCLUDED
#define mdaq_hw_server_BaseDynRegister_h_INCLUDED

#include "../HwControllerException.h"
#include <boost/cstdint.hpp>
#include <map>
#include <vector>
#include <string>

namespace mdaq
{
namespace hw
{
namespace server
{
typedef unsigned short BitOffset;
typedef unsigned short BitWidth;

/**
 * This is the base class for accessing and configuring registers dynamically.
 */
class BaseDynRegister
{
public:
	class FieldNotFound;
	class OutOfRange;

	BaseDynRegister();
	virtual ~BaseDynRegister();

	virtual void setStorageAddress(uint64_t address) = 0;
	virtual uint64_t storageAddress() const = 0;

	virtual int64_t	getValueInt64() const = 0;
	virtual void	setValueInt64(int64_t value) = 0;

	virtual int64_t	fieldReadInt64(const std::string& name) const = 0;
	virtual void	fieldWriteInt64(const std::string& name, int64_t value) = 0;

	BitOffset fieldOffset(const std::string& name) const;
	BitWidth  fieldWidth(const std::string& name) const;

	void	fieldList(std::vector<std::string>& fields) const;
	void	fieldAdd(const std::string& name, BitOffset offset, BitWidth width);
	size_t	fieldCount() const	{ return fields.size(); }
	void	fieldClear()		{ fields.clear(); }

protected:
	enum OutOfRangeType { RegisterOutOfRange, FieldOutOfRange };
	struct BitFieldParameters
	{
		BitFieldParameters(BitOffset offset, BitWidth width) : offset(offset), width(width) {}
		BitOffset offset;
		BitWidth width;
	};
	typedef std::map<std::string, BitFieldParameters> BitFieldsMap;

	BitFieldsMap fields;
};

/**
 * Exception thrown by classes deriving from BaseDynRegister when no field
 * exists for the specified name.
 */
class BaseDynRegister::FieldNotFound : public HwControllerException
{
public:
	FieldNotFound(const std::string& field);
	~FieldNotFound() throw();
	const std::string& fieldName() const { return _fieldName; }

private:
	void writeMessage(std::ostream& stream) const;
	std::string _fieldName; ///< the name of the missing field.
};

/**
 * This is a base exception class for classes deriving from BaseDynRegister to
 * signal an overflow while setting a register or a field.
 */
class BaseDynRegister::OutOfRange : public HwControllerException
{
public:
	OutOfRange(const std::string& itemName, OutOfRangeType type);
	~OutOfRange() throw() {}

protected:
	std::string itemName;
	OutOfRangeType type;
};

}
}
}

#endif // mdaq_hw_server_BaseDynRegister_h_INCLUDED
