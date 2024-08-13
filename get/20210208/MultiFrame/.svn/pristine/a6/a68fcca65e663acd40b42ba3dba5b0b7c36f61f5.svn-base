/*
 * $Id$
 * @file BitField.hpp
 * @created 28 mars 2012
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

#ifndef mfm_BitField_hpp_INCLUDED
#define mfm_BitField_hpp_INCLUDED

#include "mfm/BitField.h"
#include "mfm/Field.h"
#include <utl/BitFieldHelper.hpp>
using namespace utl;

namespace mfm {
//______________________________________________________________________
/**
 * Encodes the value of a bit field.
 * @tparam T    Type of the bit field.
 * @param val   Value to encode into the bit field.
 */
template <typename T>
void BitField::setValue(const T & val)
{
	// FIXME: users should not have to provide a type for the field containing the bit field
	T fieldValue = field_.value<T>();
	BitFieldHelper<T>::setField(fieldValue, pos_b_, width_b_, val);
	field_.setValue<T>(fieldValue);
}
//______________________________________________________________________
/**
 * Extracts the value of a bit field from its container field.
 * @tparam T    Type of the field containing the bit field.
 * @return      Value of the bit field.
 */
template <typename T>
T BitField::value() const
{
	// FIXME: users should not have to provide a type for the field containing the bit field
	T fieldValue = field_.value<T>();
	return BitFieldHelper<T>::getField(fieldValue, pos_b_, width_b_);
}
//______________________________________________________________________
} // namespace mfm
#endif /* mfm_BitField_hpp_INCLUDED */
