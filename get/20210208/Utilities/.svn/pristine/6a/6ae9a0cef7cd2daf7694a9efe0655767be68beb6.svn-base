/* =============================================================================
 * Configurator.h
 * -----------------------------------------------------------------------------
 * class utl::Configurator
 * Created on: 18 mars 2013 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef utl_Configurator_h_INCLUDED
#define utl_Configurator_h_INCLUDED

#include <utl/Singleton.hpp>
#include <utl/Exception.h>
#include <map>
#include <string>
#include <typeinfo>

namespace utl
{
/**
 * This class allows to decouple configuration aspects from a class implementation.
 * Configuration can be applied in a separate free-standing or static function that
 * will be registered in the Configurator and invoked, based on the type of the
 * object to configure and on the type of the configuration object.
 *
 * It relies on AOP (Aspect Oriented Programming) design to solve dependency
 * issues related to configuration.
 * It requires the compiler to have a working RTTI (Run-Time Type Information) implementation.
 */
class Configurator
{
	typedef void (*GenericFunction)(void*, const void*);
	typedef std::map<const char*, GenericFunction> ClassRegistry;
	typedef std::map<const char*, ClassRegistry> ConfigRegistry;
public:
	Configurator();
	~Configurator();

	template<class T, class C>
	bool registerFunction(void (*func)(T&, const C&));

	template<class T, class C>
	Configurator& add(void (*func)(T&, const C&));

	template<class T, class C>
	void configure(T& object, const C& config);

	class MissingEntry : public Exception
	{
	public:
		enum Type {Config, Class};

		MissingEntry(Type t, const char* typeName);

		void writeMessage(std::ostream& stream) const;

		const Type type;
		const char* const typeName;
	};

private:
	bool setFunction(const char* configTypeId, const char* objTypeId, GenericFunction function);
	GenericFunction getFunction(const char* configTypeId, const char* objTypeId);

	ConfigRegistry registry;
};

typedef ::utl::Singleton<Configurator> GlobalConfigurator;

#define REGISTER_CONFIGURATOR(function) \
	static bool BOOST_JOIN(configuratorRegistered, __LINE__) = ::utl::GlobalConfigurator::instance().registerFunction(function);

// =============================================================================
// utl::Configurator inline class methods body
// -----------------------------------------------------------------------------

/**
 * Register a configuration function.
 * It must be a static or free-standing function taking 2 arguments: the first
 * is a reference to the object to configure, the second is a const reference to
 * the configuration to apply.
 */
template<class T, class C>
bool Configurator::registerFunction(void (*func)(T&, const C&))
{
	return setFunction(typeid(C).name(), typeid(T).name(), reinterpret_cast<GenericFunction>(func));
}

template<class T, class C>
Configurator& Configurator::add(void (*func)(T&, const C&))
{
	registerFunction(func);
	return *this;
}

/**
 * Configures an object using the previously registered function matching the type
 * of the object and the type of the configuration.
 */
template<class T, class C>
void Configurator::configure(T& object, const C& config)
{
	typedef void (*Function)(T&, const C&);
	Function configurator =
		reinterpret_cast<Function>(
			getFunction(typeid(config).name(), typeid(object).name())
		);
	configurator(object, config);
}

} /* namespace utl */

#endif /* utl_Configurator_h_INCLUDED */
