/**
 * @file DataLinkTest.cpp
 * @date 11 mars 2013
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetBench-GANIL software project.
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

#include "get/rc/DataLinkSet.h"
#include <boost/test/auto_unit_test.hpp>
#include <iostream>

namespace get
{
namespace rc
{

BOOST_AUTO_TEST_CASE(dataSenderToXml_test)
{
	DataSenderId sender("CoBo[Crate0Slot0].Asad[0]");

	const std::string expectedString = "<DataSender id=\"CoBo[Crate0Slot0].Asad[0]\" />";

	const std::string actualString = sender.toXml();

	BOOST_CHECK_EQUAL(actualString, expectedString);
}

BOOST_AUTO_TEST_CASE(xmlToDataSender_test)
{
	DataSenderId sender;

	const std::string initialString = "<DataSender id=\"Mutant\"></DataSender>";

	BOOST_CHECK_EQUAL(sender.fromXml(initialString), true);
	BOOST_CHECK_EQUAL(sender.id, "Mutant");
}

BOOST_AUTO_TEST_CASE(xmlToDataLink_test)
{
	const std::string initialString = "<DataLink>\n"
			"\t<DataSender id=\"CoBo[Crate0Slot1].Asad[3]\"/>\n"
			"\t<DataRouter name=\"Narval3\" ipAddress=\"132.166.9.3\" port=\"46005\" type=\"FDT\"/>\n"
			"</DataLink>";

	DataLink link;

	BOOST_CHECK_EQUAL(link.fromXml(initialString), true);
	BOOST_CHECK_EQUAL(link.sender.id, "CoBo[Crate0Slot1].Asad[3]");
	BOOST_CHECK_EQUAL(link.router.name, "Narval3");
	BOOST_CHECK_EQUAL(link.router.sockAddr.ipAddress().toString(), "132.166.9.3");
	BOOST_CHECK_EQUAL(link.router.sockAddr.port(), 46005u);
	BOOST_CHECK_EQUAL(link.router.flowType, "FDT");
}

BOOST_AUTO_TEST_CASE(xmlToDataLinkSet_test)
{
	const std::string initialString = "<DataLinkSet><DataLink><DataSender id=\"CoBo[Crate00_Slot00].AsAd[0]\""
			" /><DataRouter ipAddress=\"132.166.9.203\" port=\"46005\" type=\"Tcp\" /></DataLink></DataLinkSet>";

	DataLinkSet links;

	BOOST_CHECK_EQUAL(links.fromXml(initialString), true);
	BOOST_CHECK_EQUAL(links.size(), 1u);
	DataSenderId sender("CoBo[Crate00_Slot00].AsAd[0]");
	BOOST_CHECK_EQUAL((*links.begin()).sender == sender, true);
}

} /* namespace rc */
} /* namespace get */


