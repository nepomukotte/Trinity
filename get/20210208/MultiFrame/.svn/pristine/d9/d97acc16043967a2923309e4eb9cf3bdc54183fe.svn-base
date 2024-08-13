/*
 * $Id$
 * @file Field.hpp
 * @created 20 févr. 2012
 * @author sizun
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 * -----------------------------------------------------------------------------
 */

#ifndef mfm_Field_hpp_INCLUDED
#define mfm_Field_hpp_INCLUDED

#include <mfm/Field.h>
#include <mfm/AbstractFieldContainer.h>
#include <utl/BinIO.hpp>
#include <utl/BitFieldHelper.hpp>
using namespace utl;

namespace mfm {
//______________________________________________________________________
/**
 * Encodes the value of a word of type T into the size of the field and according to its endianness.
 * @tparam T Type of the input word.
 * @param value Input word from which bytes are encoded.
 */
template <typename T>
void Field::setValue(const T & value)
{
	BinIO::encode<T>(container_->header().endianness(), serializer_.size_B(), value, serializer_.outputStream());
}
//______________________________________________________________________
/**
 * Decodes the field according to its endianness into a word of type T.
 * @tparam T Type of the input word.
 * @return Decoded value of the field.
 */
template <typename T>
T Field::value() const
{
	T value;
	BinIO::decode<T>(container_->header().endianness(), serializer_.size_B(), value, serializer_.inputStream());
	return value;
}
//______________________________________________________________________
} // namespace mfm
#endif // mfm_Field_hpp_INCLUDED
