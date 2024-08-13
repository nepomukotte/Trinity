/* =============================================================================
 * UdpTransceiverTest.cpp
 * -----------------------------------------------------------------------------
 * class ::
 * Created on: 10 déc. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This software is part of the MDaq project.
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
 * of their systems and/or ref to be ensured and, more generally, to use and
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
#include <mdaq/utl/net/UdpTransceiver.h>
#include <mdaq/mt/SysClock.h>
#include <utl/Logging.h>
#include <boost/test/auto_unit_test.hpp>

using namespace utl::net;
using namespace mdaq::utl::net;

struct UdpTransceiverFixture : NetTestClient
{
	UdpTransceiverFixture()
		: NetTestClient(), sender(),
		  receiver(SocketAddress(IpAddress(127, 0, 0, 1), 10003)),
		  serverAddr(IpAddress(127, 0, 0, 1), 10002)
	{
		mdaq::mt::SysClock::sleep_us(100000);
	}
	~UdpTransceiverFixture()
	{
	}

	NetTestData ref;
	UdpTransceiver sender;
	UdpTransceiver receiver;
	SocketAddress serverAddr;
};

BOOST_FIXTURE_TEST_SUITE(mdaq_utl_net_UdpTransceiverTest, UdpTransceiverFixture)

BOOST_AUTO_TEST_CASE(udp_send1)
{
	sender.doSend(serverAddr, buffer(ref.dataset1));
}

BOOST_AUTO_TEST_CASE(udp_send2)
{
	BOOST_CHECK_THROW(sender.doSend(serverAddr, buffer(ref.dataset2)), UdpTransceiver::SendException);
	sender.doSend(serverAddr, buffer(ref.dataset1));
}

BOOST_AUTO_TEST_CASE(udp_recv1)
{
	NetTestData::dataset3_type data(ref.dataset3.size());
	SocketAddress source;
	receiver.doRcve(source, buffer(data), 500);
	BOOST_CHECK(data == ref.dataset3);
	BOOST_CHECK_EQUAL(source.ipAddress(), serverAddr.ipAddress());
}

BOOST_AUTO_TEST_CASE(udp_recv_timeout)
{
	NetTestData::dataset1_type data(ref.dataset1.size());
	SocketAddress source;
	BOOST_CHECK_THROW(receiver.doRcve(source, buffer(data), 500), Socket::TimeoutException);
}

BOOST_AUTO_TEST_SUITE_END()
