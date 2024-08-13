/* =============================================================================
 * TcpConnectionTest.cpp
 * -----------------------------------------------------------------------------
 * class ::
 * Created on: 4 déc. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include "NetTestClient.h"
#include "NetTestData.h"
#include <utl/Logging.h>
#include <boost/test/auto_unit_test.hpp>

using namespace utl::net;
using namespace mdaq::utl::net;

struct TcpConnectionFixture : NetTestClient
{
	TcpConnectionFixture()
		: NetTestClient()
	{
	}
	~TcpConnectionFixture()
	{
	}

	NetTestData data;
};

BOOST_FIXTURE_TEST_SUITE(mdaq_utl_net_TcpConnectionTest, TcpConnectionFixture)

BOOST_AUTO_TEST_CASE(tcp_send1)
{
	conn.doSend(buffer(data.dataset1), SEND_DATA_TIMEOUT);
}

BOOST_AUTO_TEST_CASE(tcp_send2)
{
	conn.doSend(buffer(data.dataset2), SEND_DATA_TIMEOUT);
}

BOOST_AUTO_TEST_CASE(tcp_send4)
{
	conn.doSend(buffer(data.dataset4), SEND_DATA_TIMEOUT);
}

BOOST_AUTO_TEST_CASE(tcp_send5)
{
	conn.doSend(buffer(data.dataset5), SEND_DATA_TIMEOUT);
}

BOOST_AUTO_TEST_CASE(tcp_send6)
{
	conn.doSend(buffer(data.dataset6), SEND_DATA_TIMEOUT);
}

//BOOST_AUTO_TEST_CASE(test_send_timeout)
//{
//	BOOST_CHECK_THROW(conn.doSend(&data.dataset2[0], data.dataset2.size() * sizeof(data.dataset2[0]), 500), Socket::TimeoutException);
//}

BOOST_AUTO_TEST_SUITE_END()
