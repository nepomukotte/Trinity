/**
 * @file EnumStringMapper.cpp
 * @date 21 juil. 2014
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

#include "EnumStringMapper.hpp"

namespace utl
{

// =============================================================================
// utl::EnumStringNotFound methods body
// -----------------------------------------------------------------------------

EnumStringNotFound::EnumStringNotFound(const std::type_info& eti, const std::string& valueStr)
	: valueStr(valueStr), enumName(eti.name())
{
}

EnumStringNotFound::~EnumStringNotFound() throw()
{
}

void EnumStringNotFound::writeMessage(std::ostream& stream) const
{
	stream << "Cannot find an enumerated value named '" << valueStr << "' in enum " << enumName;
}

// =============================================================================
// utl::EnumMapperBadDefinition methods body
// -----------------------------------------------------------------------------

EnumMapperBadDefinition::EnumMapperBadDefinition(const std::type_info& eti)
{
	msgBuilder << "Bad definition of enum: " << eti.name() << ' ';
}

EnumMapperBadDefinition::EnumMapperBadDefinition(const EnumMapperBadDefinition& rvalue)
	: Exception(rvalue)
{
	msgBuilder.str(rvalue.msgBuilder.str());
}

EnumMapperBadDefinition::~EnumMapperBadDefinition() throw()
{
}

EnumMapperBadDefinition& EnumMapperBadDefinition::operator=(const EnumMapperBadDefinition& rvalue)
{
	Exception::operator=(rvalue);
	msgBuilder.str(rvalue.msgBuilder.str());
	return *this;
}

void EnumMapperBadDefinition::writeMessage(std::ostream& stream) const
{
	stream << msgBuilder.str();
}

} /* namespace utl */
