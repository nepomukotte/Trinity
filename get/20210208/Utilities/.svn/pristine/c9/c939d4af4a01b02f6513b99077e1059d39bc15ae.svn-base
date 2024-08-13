/**
 * @file Construct.hpp
 * @date 1 août 2013
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

#ifndef utl_Construct_hpp_INCLUDED
#define utl_Construct_hpp_INCLUDED

#include <boost/preprocessor/inc.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/repetition/enum_shifted.hpp>
#include <boost/preprocessor/repetition/enum_shifted_params.hpp>
#include <boost/pointee.hpp>
#include <boost/type_traits/remove_cv.hpp>

#define CONSTRUCT_MAX_ARITY 3

#define DEC_ARITY(r, arity)     BOOST_PP_DEC(arity)
#define INC_ARITY(r, arity)     BOOST_PP_INC(arity)
#define NOT_ZERO(r, i)          BOOST_PP_NOT_EQUAL(i, 0)
#define NOT_MAX_ARITY(r, state) BOOST_PP_NOT_EQUAL(state, BOOST_PP_INC(CONSTRUCT_MAX_ARITY))

#define ARG_TYPE(i)             BOOST_PP_CAT(T,i)

#define TEMPLATE_ITEM(J,I,D)    class ARG_TYPE(I)
#define TEMPLATE_LIST(arity)    BOOST_PP_ENUM_SHIFTED(BOOST_PP_INC(arity),TEMPLATE_ITEM,BOOST_PP_EMPTY)

#define FUNCTION_PARM(J,I,D)    const ARG_TYPE(I)& BOOST_PP_CAT(arg,I)
#define FUNCTION_PARMS(arity)   BOOST_PP_ENUM_SHIFTED(BOOST_PP_INC(arity),FUNCTION_PARM,BOOST_PP_EMPTY)
#define FUNCTION_ARGS(arity)    BOOST_PP_ENUM_SHIFTED_PARAMS(BOOST_PP_INC(arity), arg)

#define CALLER_IMPL(r, arity) \
	template<TEMPLATE_LIST(arity)> \
	result_type operator()(FUNCTION_PARMS(arity)) const { return make_pointer(new value_type(FUNCTION_ARGS(arity))); } \
	\

namespace utl
{
template<class Pointer>
struct Construct
{
        typedef typename boost::remove_cv<Pointer>::type result_type;
        typedef typename boost::pointee<result_type>::type value_type;

		result_type make_pointer(value_type* rawptr) const
		{
			return result_type(rawptr);
		}

        result_type operator()() const { return make_pointer(new value_type()); }
		BOOST_PP_FOR(1, NOT_MAX_ARITY, INC_ARITY, CALLER_IMPL)
};
} /* namespace utl */

#undef CONSTRUCT_MAX_ARITY

#undef DEC_ARITY
#undef INC_ARITY
#undef NOT_ZERO
#undef NOT_MAX_ARITY
#undef ARG_TYPE
#undef TEMPLATE_ITEM
#undef TEMPLATE_LIST
#undef FUNCTION_PARM
#undef FUNCTION_PARMS
#undef FUNCTION_ARGS
#undef CALLER_IMPL

#endif /* utl_Construct_hpp_INCLUDED */
