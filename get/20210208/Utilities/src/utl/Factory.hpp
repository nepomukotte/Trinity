/* =============================================================================
 * Factory.hpp
 * -----------------------------------------------------------------------------
 * class utl::Factory
 * class utl::FactoryMap
 * class utl::MissingFactory
 * Created on: 27 nov. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef utl_Factory_hpp_INCLUDED
#define utl_Factory_hpp_INCLUDED

#include <utl/Exception.h>
#include <boost/function.hpp>
#include <boost/pointee.hpp>
#include <boost/static_assert.hpp>
#include <boost/preprocessor/inc.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/repetition/enum_shifted.hpp>
#include <boost/preprocessor/repetition/enum_shifted_params.hpp>
#include <map>
#include <memory>

/*
 * Redefine FACTORY_MAX_ARITY, prior to include this header,
 * if you want to pass more parameters to your factories.
 */
#ifndef FACTORY_MAX_ARITY
#define FACTORY_MAX_ARITY 3
#endif

#define DEC_ARITY(r, arity)     BOOST_PP_DEC(arity)
#define INC_ARITY(r, arity)     BOOST_PP_INC(arity)
#define NOT_ZERO(r, i)          BOOST_PP_NOT_EQUAL(i, 0)
#define NOT_MAX_ARITY(r, state) BOOST_PP_NOT_EQUAL(state, BOOST_PP_INC(FACTORY_MAX_ARITY))

#define ARG_TYPE(i)             BOOST_PP_CAT(BOOST_PP_CAT(arg,i),_type)
#define ARG_TYPEDEF(r, state)   typedef typename Base::FactoryMethod::ARG_TYPE(state) ARG_TYPE(state);

#define FUNCTION_PARM(J,I,D)    ARG_TYPE(I) BOOST_PP_CAT(arg,I)
#define FUNCTION_PARMS(arity)   BOOST_PP_ENUM_SHIFTED(BOOST_PP_INC(arity),FUNCTION_PARM,BOOST_PP_EMPTY)
#define FUNCTION_ARGS(arity)    BOOST_PP_ENUM_SHIFTED_PARAMS(BOOST_PP_INC(arity), arg)

#define FACTORY_IMPL_SPEC(r, arity) \
	template<typename Base> \
	struct FactoryImpl<Base, arity> : Base \
	{ \
		typedef typename Base::ProductPtr ProductPtr; \
		BOOST_PP_FOR(arity, NOT_ZERO, DEC_ARITY, ARG_TYPEDEF) \
		ProductPtr create(FUNCTION_PARMS(arity)) const { return ProductPtr(this->method_(FUNCTION_ARGS(arity))); } \
	};

namespace utl
{
/**
 * This is a generic class whose specializations allow to define a create()
 * method for any possible number of parameters (arity).
 */
template<typename Base, unsigned arity>
struct FactoryImpl
{
	BOOST_STATIC_ASSERT(arity <= FACTORY_MAX_ARITY);
};

/*
 * Generating specializations for arity (ie: number of parameters) between 0 and FACTORY_MAX_ARITY.
 */
BOOST_PP_FOR(0, NOT_MAX_ARITY, INC_ARITY, FACTORY_IMPL_SPEC)

template<typename Signature>
class BaseFactory
{
public:
	typedef boost::function<Signature> FactoryMethod;
	typedef typename FactoryMethod::result_type ReturnType;
	typedef typename boost::pointee<ReturnType>::type ProductType;
	typedef std::auto_ptr<ProductType> ProductPtr;

	BOOST_STATIC_ASSERT(boost::is_pointer<ReturnType>::value);

	void setMethod(const FactoryMethod& newFactory) { method_ = newFactory; }
	const FactoryMethod& method() const { return method_; }

protected:
	FactoryMethod method_;
};

/**
 * This class allows to implement abstract factory pattern easily.
 * @tparam Signature the signature of the factory function. Currently, it must
 * return a raw pointer. Ownership of this pointer will be taken by an auto_ptr
 * which will be returned by create(). The factory function can have any
 * number of parameters, the create() method will be generated with the same
 * parameters.
 *
 * Compared to other possible implementations:
 *  - unlike traditional GOF pattern it doesn't require you to define classes
 *    for abstract and concrete factories: a method is enough (and with
 *    boost::factory nothing is needed)
 *  - unlike function pointer solution it allows you to use any kind of callable:
 *    free standing function, method, functor, using argument binding if necessary (boost::bind)
 *  - unlike using boost::function directly, it handles ownership and smart-pointers
 *    upcasts without problem, and have an interface that makes explicit the
 *    usage as a factory.
 */
template<typename Signature>
struct Factory : FactoryImpl<BaseFactory<Signature>, boost::function<Signature>::arity>
{
};

/**
 * This class manages factories all creating products of the same abstract
 * type, each factory being registered and associated to a key.
 */
template<class Key, typename Signature>
class FactoryMap
{
public:
	typedef Key FactoryKey;
	typedef typename BaseFactory<Signature>::FactoryMethod FactoryMethod;

	FactoryMap<Key, Signature>& add(const Key& key, const FactoryMethod& factory);
	const Factory<Signature>& operator[](const Key& key) const;

	bool hasFactory(const Key& key);

protected:
	typedef std::map<Key, Factory<Signature> > Map;
	Map factories;
};

/**
 * This exception reports the FactoryMap has been unable to find the requested
 * factory.
 */
class MissingFactory : public Exception
{
public:
	template<class Key>
	MissingFactory(const Key& key, const std::type_info& type);
	~MissingFactory() throw();

	void writeMessage(std::ostream& stream) const;

	std::string factoryKey;
	std::string abstractType;
};

/**
 * This class allows to register a factory by creating a global static instance
 * of it. This is further simplified by the REGISTER_FACTORY macro.
 */
template<class F>
class FactoryRegistrar
{
public:
	typedef typename F::elem_type               FactoryType;
	typedef typename FactoryType::FactoryKey    FactoryKey;
	typedef typename FactoryType::FactoryMethod FactoryMethod;

	FactoryRegistrar(const FactoryKey& key, const FactoryMethod& method);
	~FactoryRegistrar();
};

} /* namespace utl */

#define REGISTER_FACTORY(factory, key, method) \
	static ::utl::FactoryRegistrar<factory> BOOST_JOIN(fr, __LINE__)(key, method);

#undef FACTORY_IMPL_SPEC
#undef FUNCTION_ARGS
#undef FUNCTION_PARMS
#undef FUNCTION_PARM
#undef ARG_TYPEDEF
#undef ARG_TYPE
#undef NOT_MAX_ARITY
#undef NOT_ZERO
#undef INC_ARITY
#undef DEC_ARITY
#undef FACTORY_MAX_ARITY

#endif /* utl_Factory_hpp_INCLUDED */
