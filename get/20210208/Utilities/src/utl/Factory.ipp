/* =============================================================================
 * Factory.cpp
 * -----------------------------------------------------------------------------
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

#ifndef utl_Factory_ipp_INCLUDED
#define utl_Factory_ipp_INCLUDED

#include "Factory.hpp"
#include <boost/lexical_cast.hpp>

namespace utl
{
// =============================================================================
// utl::FactoryMap class methods body
// -----------------------------------------------------------------------------

template<class Key, typename Signature>
FactoryMap<Key, Signature>& FactoryMap<Key, Signature>::add(const Key& key, const FactoryMethod& factory)
{
	factories[key].setMethod(factory);
	return *this;
}

template<class Key, typename Signature>
const Factory<Signature>& FactoryMap<Key, Signature>::operator[](const Key& key) const
{
	typename Map::const_iterator f = factories.find(key);
	if(f != factories.end())
	{
		return f->second;
	}
	else
	{
		throw MissingFactory(key, typeid(typename BaseFactory<Signature>::ProductType));
	}
}

template<class Key, typename Signature>
bool FactoryMap<Key, Signature>::hasFactory(const Key& key)
{
	typename Map::const_iterator f = factories.find(key);
	return f != factories.end();
}

// =============================================================================
// utl::MissingFactory class methods body
// -----------------------------------------------------------------------------

template<class Key>
MissingFactory::MissingFactory(const Key& key, const std::type_info& type)
	: factoryKey(boost::lexical_cast<std::string>(key)),
	  abstractType(type.name())
{
}

inline void MissingFactory::writeMessage(std::ostream& stream) const
{
	stream << "Cannot find a factory for type \'" << abstractType << "\' of kind \'" << factoryKey << '\'';
}

inline MissingFactory::~MissingFactory() throw()
{
}

// =============================================================================
// utl::FactoryRegistrar class methods body
// -----------------------------------------------------------------------------

template<class F>
FactoryRegistrar<F>::FactoryRegistrar(const FactoryKey& key, const FactoryMethod& method)
{
	F::instance().add(key, method);
}

template<class F>
FactoryRegistrar<F>::~FactoryRegistrar()
{
}

} /* namespace utl */

#endif /* utl_Factory_ipp_INCLUDED */
