/**
 * @file ArrayHelper.hpp
 * @date 27 août 2013
 * @author fchateau
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Frederic Chateau
 *
 * This file is part of the Utilities software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#ifndef utl_ArrayHelper_hpp_INCLUDED
#define utl_ArrayHelper_hpp_INCLUDED

#include <algorithm>

namespace utl
{
namespace array
{
   template<typename T, int N>
   size_t size(T (&)[N])
   {
	  return N;
   }
   template<typename T, int N>
   T* begin(T (&a)[N])
   {
	  return &a[0];
   }
   template<typename T, int N>
   T* end(T (&a)[N])
   {
	  return &a[N];
   }

   template<class InputIterator, class T>
   size_t find_index(InputIterator first, InputIterator last, const T& val)
   {
	   return std::distance(first, std::find(first, last, val));
   }

} /* namespace array */
} /* namespace utl */

#endif /* utl_ArrayHelper_hpp_INCLUDED */
