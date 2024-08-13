/**
 * @file algorithm.h
 * @date Nov 18, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: algorithm.h 1015 2013-11-18 15:06:20Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the CoBoFrameViewer software project.
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

#ifndef get_ana_algorithm_h_INCLUDED
#define get_ana_algorithm_h_INCLUDED

#include "utl/Logging.h"
#include <algorithm>
#include <cstddef>

namespace get {
namespace ana {
namespace algorithm {
//_________________________________________________________________________________________________
/**
 * Computes peaking time of a data sample.
 */
template < class ForwardIter, typename T >
size_t peakingTime(ForwardIter first, ForwardIter last, const T & threshold, const T & offset)
{
	// Find bucket and value of max sample
	ForwardIter maxIter = std::max_element(first, last);
	const T deltaMax = *maxIter - offset;

	// Find first value, left of max, greater than threshold*max
	ForwardIter leftIter = std::find_if(first, maxIter, std::bind2nd(std::greater< uint16_t >(), offset + threshold*deltaMax));

	//LOG_WARN() << "5pc=" << *leftIter <<  " max=" << *maxIter;

	// Compute full width
	return std::distance(leftIter, maxIter);
}
//_________________________________________________________________________________________________
/**
 * Computes FWHM (Full Width at Half-Maximum) of a data sample.
 */
template < class ForwardIter, typename T >
size_t fwhm(ForwardIter first, ForwardIter last, const T & threshold, const T & offset)
{
	// Find bucket and value of max sample
	ForwardIter maxIter = std::max_element(first, last);
	const T deltaMax = *maxIter - offset;

	// Find first value, left of max, greater than half-max
	ForwardIter leftIter = std::find_if(first, maxIter, std::bind2nd(std::greater_equal< uint16_t >(), offset + threshold*deltaMax));

	// Find first value, right of max, less than half-max
	ForwardIter rightIter = std::find_if(maxIter+1, last, std::bind2nd(std::less< uint16_t >(), offset + threshold*deltaMax));

	//LOG_WARN() << "left=" << *leftIter <<  " max=" << *maxIter << " right=" << *rightIter;

	// Compute full width
	return std::distance(leftIter, rightIter);
}
//_________________________________________________________________________________________________
} // namespace algorithm
} // namespace ana
} //namespace get

#endif /* get_ana_algorithm_h_INCLUDED */
