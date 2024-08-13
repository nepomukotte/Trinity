/* =============================================================================
 * Configurator.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::utl::Configurator
 * Created on: 18 mars 2013 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This software is part of the MDaq2 project.
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

#include "Configurator.h"
#include <utl/Singleton.ipp>
#include <stdexcept>

namespace utl
{

Configurator::Configurator()
{
}

Configurator::~Configurator()
{
}

bool Configurator::setFunction(const char* configTypeId, const char* objTypeId, GenericFunction function)
{
	ClassRegistry& classRegistry = registry[configTypeId];
	std::pair<ClassRegistry::iterator, bool> result = classRegistry.insert(std::make_pair(objTypeId, function));
	return result.second;
}

Configurator::GenericFunction Configurator::getFunction(const char* configTypeId, const char* objTypeId)
{
	ConfigRegistry::iterator itConfig = registry.find(configTypeId);
	if(itConfig == registry.end())
	{
		throw MissingEntry(MissingEntry::Config, configTypeId);
	}
	ClassRegistry& classRegistry = itConfig->second;
	ClassRegistry::iterator itClass = classRegistry.find(objTypeId);
	if(itClass == classRegistry.end())
	{
		throw MissingEntry(MissingEntry::Class, objTypeId);
	}
	return itClass->second;
}

template class Singleton<Configurator>;

// =============================================================================
// utl::Configurator::MissingEntry class methods body
// -----------------------------------------------------------------------------

Configurator::MissingEntry::MissingEntry(Type t, const char* name)
	: type(t), typeName(name)
{
}

void Configurator::MissingEntry::writeMessage(std::ostream& stream) const
{
	stream
		<< "Missing "
		<< (type == Config ? "configuration" : "class") << " entry '"
		<< typeName << "' in configurator.";
}

} /* namespace utl */
