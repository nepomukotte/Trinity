/**
 * @file ConfigIdSetTest.cpp
 * @date 15 févr. 2013
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetBench software project.
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


#include "get/rc/ConfigIdSet.h"
#include <boost/test/auto_unit_test.hpp>
#include <iostream>

namespace get
{
namespace rc
{

BOOST_AUTO_TEST_CASE(configIdSetToXml_test)
{
	ConfigIdSet configIdSet;
	configIdSet.insert("monoNode", "prepareParam", "pedestals");

	const std::string expectedString =
			"<ConfigIdSet>\n"
			"	<ConfigId>"
				"<SubConfigId type=\"describe\">monoNode</SubConfigId>"
				"<SubConfigId type=\"prepare\">prepareParam</SubConfigId>"
				"<SubConfigId type=\"configure\">pedestals</SubConfigId>"
			"</ConfigId>\n"
			"</ConfigIdSet>";

	const std::string actualString = configIdSet.toXml();

	BOOST_CHECK_EQUAL(actualString, expectedString);
}

BOOST_AUTO_TEST_CASE(xmlToConfigIdSet_test)
{
	const std::string initialString =
			"<ConfigIdSet>\n"
			"<ConfigId>\n"
				"\t<SubConfigId type=\"describe\">monoNode</SubConfigId>\n"
				"\t<SubConfigId type=\"prepare\">prepareParam</SubConfigId>\n"
				"\t<SubConfigId type=\"configure\">pedestals</SubConfigId>\n"
			"</ConfigId>\n"
			"</ConfigIdSet>";

	ConfigIdSet configIdSet;

	BOOST_CHECK_EQUAL(configIdSet.fromXml(initialString), true);

	BOOST_CHECK_EQUAL(configIdSet.size(), 1u);

	ConfigId configId("monoNode", "prepareParam", "pedestals");
	BOOST_CHECK_EQUAL(configIdSet.find(configId) == configIdSet.end(), false);
}

} /* namespace rc */
} /* namespace get */



