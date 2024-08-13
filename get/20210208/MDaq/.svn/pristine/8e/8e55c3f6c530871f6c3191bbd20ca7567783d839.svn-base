/* =============================================================================
 * BaseDynRegister.cpp
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

#include "BaseDynRegister.h"

namespace mdaq
{
namespace hw
{
namespace server
{
/**
 * Constructor.
 */
BaseDynRegister::BaseDynRegister()
{
}

/**
 * Destructor.
 */
BaseDynRegister::~BaseDynRegister()
{
}

/**
 * @fn int64_t BaseDynRegister::getValueInt64() const
 * Returns the value of the register whatever its actual size provided it does not exceed 64 bits.
 */

/**
 * @fn void BaseDynRegister::setValueInt64(int64_t value)
 * Sets the value of the register whatever its actual size provided it does not exceed 64 bits.
 */

/**
 * @fn int64_t BaseDynRegister::fieldReadInt64(const std::string& name) const
 * Returns the value of the field having the specified name.
 */

/**
 * @fn void BaseDynRegister::fieldWriteInt64(const std::string& name, int64_t value)
 * Sets the value of the field having the specified name.
 */

/**
 * @fn size_t BaseDynRegister::fieldCount() const
 * Returns the number of fields composing the register.
 */

/**
 * @fn void BaseDynRegister::fieldClear()
 * Clears the list of fields composing the register.
 */

/**
 * Returns offset of specified field.
 * @param name Name of the field.
 * @return Field offset.
 */
BitOffset BaseDynRegister::fieldOffset(const std::string& name) const
{
	BitFieldsMap::const_iterator it = fields.find(name);
	if (it == fields.end())
		throw FieldNotFound(name);
	const BaseDynRegister::BitFieldParameters& field = it->second;
	return field.offset;
}

/**
 * Returns width of specified field.
 * @param name Name of the field.
 * @return Field width.
 */
BitWidth  BaseDynRegister::fieldWidth(const std::string& name) const
{
	BitFieldsMap::const_iterator it = fields.find(name);
	if (it == fields.end())
		throw FieldNotFound(name);
	const BaseDynRegister::BitFieldParameters& field = it->second;
	return field.width;
}

/**
 * Fills a string vector with the list of the names of the register's fields.
 */
void BaseDynRegister::fieldList(std::vector<std::string>& fieldList) const
{
	for(BitFieldsMap::const_iterator it = fields.begin(), end = fields.end();
		it != end; ++it)
	{
		fieldList.push_back(it->first);
	}
}

/**
 * Adds a new bit-field to this register.
 * @param name the name of the field. It will be used to reference it.
 * @param offset the position of the first bit.
 * @param width the number of bit composing the field.
 */
void BaseDynRegister::fieldAdd(const std::string& name, BitOffset offset, BitWidth width)
{
	fields.insert(std::make_pair(name, BitFieldParameters(offset, width)));
}


//____________________________________________________________________________//

/**
 * Exception constructor.
 * @param fieldName The name of the missing field.
 */
BaseDynRegister::FieldNotFound::FieldNotFound(const std::string& fieldName)
	: _fieldName(fieldName)
{
}

/**
 * Destructor.
 */
BaseDynRegister::FieldNotFound::~FieldNotFound() throw()
{
}

/**
 * Implementation of formatted message as expected by base class std::exception.
 */
void BaseDynRegister::FieldNotFound::writeMessage(std::ostream& stream) const
{
	stream << "No field is named \"" << fieldName() << '"';
}

//____________________________________________________________________________//

/**
 * Exception constructor.
 */
BaseDynRegister::OutOfRange::OutOfRange(const std::string& itemName, OutOfRangeType type)
	: itemName(itemName),
	type(type)
{
}

}
}
}
