/* =============================================================================
 * ForwardCtor.hpp
 * -----------------------------------------------------------------------------
 * class ::
 * Created on: 7 déc. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This software is part of the Utilities project.
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

#ifndef utl_ForwardCtor_hpp_INCLUDED
#define utl_ForwardCtor_hpp_INCLUDED

#include <boost/utility/enable_if.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/repetition/enum_shifted_params.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/inc.hpp>
#include <boost/version.hpp>

#ifndef UTL_FORWARD_CTOR_MAX_ARITY
#define UTL_FORWARD_CTOR_MAX_ARITY 5
#endif

#define INC_ARITY(r, arity)     BOOST_PP_INC(arity)
#define NOT_MAX_ARITY(r, state) BOOST_PP_NOT_EQUAL(state, BOOST_PP_INC(UTL_FORWARD_CTOR_MAX_ARITY))

#define FORWARD_CTOR_IMPL(r,N) \
	template<BOOST_PP_ENUM_PARAMS(N,class A)> \
	ForwardCtor(const BOOST_PP_CAT(boost::in_place_factory,N)<BOOST_PP_ENUM_PARAMS(N,A)>& f) : T(BOOST_PP_ENUM_PARAMS(N,f.m_a)) {}

namespace utl
{
/**
 * This class can forward any number of parameters to a base-class constructor.
 * It is based on boost::in_place factory: the derived class must inherit from
 * ForwardCtor<Base>, "Base" being the class to which we want to forward
 * parameters.
 * The derived class also must declare a constructor with a template parameter:
 * the in_place factory wrapping the parameters to forward.
 */
template<class T>
struct ForwardCtor : T
{
	template<class Factory>
	ForwardCtor(const Factory&, typename boost::enable_if<boost::is_same<Factory, boost::in_place_factory0>, int>::type = 0) : T() {}

	BOOST_PP_FOR(1, NOT_MAX_ARITY, INC_ARITY, FORWARD_CTOR_IMPL)
};

} /* namespace utl */

#undef INC_ARITY
#undef NOT_MAX_ARITY
#undef FORWARD_CTOR_IMPL
#undef UTL_FORWARD_CTOR_MAX_ARITY

#endif /* utl_ForwardCtor_hpp_INCLUDED */
