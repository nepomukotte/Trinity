/* =============================================================================
 * TestServer.cpp
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

#include "NetTestData.h"
#include <mdaq/utl/net/TcpAcceptor.h>
#include <mdaq/utl/net/UdpTransceiver.h>
#include <mdaq/mt/SysClock.h>
#include <utl/Logging.h>
#include <map>

using namespace utl::net;
using namespace mdaq::utl::net;

static const Milliseconds ACCEPT_TIMEOUT = INDEFINITE;
static const Milliseconds SEND_REPLY_TIMEOUT = 500;
static const Milliseconds SEND_DATA_TIMEOUT = 2000;
static const Milliseconds RCVE_DATA_TIMEOUT = 2000;

template<class T, size_t N>
bool arrayEquality(const T (&ref)[N], const T (&data)[N])
{
	bool isEqual = true;
	for(size_t i = 0; i < N and isEqual; ++i)
	{
		if(ref[i] != data[i]) isEqual = false;
	}
	return isEqual;
}

template<class T>
void checkResult(TcpConnection& conn, const T& ref, const T& other)
{
	std::string result = (ref == other) ? "Ok" : "Failed";
	conn.doSendValue(result, SEND_REPLY_TIMEOUT);
}

template<class T>
void tcpSendTest(TcpConnection& conn, const T& ref)
{
	T data(ref.size(), 0);
	conn.doRcve(buffer(data), RCVE_DATA_TIMEOUT);
	checkResult(conn, ref, data);
}

void tcp_send1(TcpConnection& conn)
{
	NetTestData ref;
	tcpSendTest(conn, ref.dataset1);
}

void tcp_send2(TcpConnection& conn)
{
	NetTestData ref;
	tcpSendTest(conn, ref.dataset2);
}

void tcp_send4(TcpConnection& conn)
{
	NetTestData ref;
	tcpSendTest(conn, ref.dataset4);
}

void tcp_send5(TcpConnection& conn)
{
	NetTestData ref;
	NetTestData::dataset5_type data = {0};
	conn.doRcve(buffer(data), RCVE_DATA_TIMEOUT);
	std::string result = arrayEquality(ref.dataset5, data) ? "Ok" : "Failed";
	conn.doSendValue(result, SEND_REPLY_TIMEOUT);
}

void tcp_send6(TcpConnection& conn)
{
	NetTestData ref;
	NetTestData::dataset6_type data;
	conn.doRcve(buffer(data), RCVE_DATA_TIMEOUT);
	std::string result = (ref.dataset6 == data) ? "Ok" : "Failed";
	conn.doSendValue(result, SEND_REPLY_TIMEOUT);
}

void tcp_send_timeout(TcpConnection& conn)
{
	//usleep(750000);
	//test_send2(conn);
	conn.doSendValue(std::string("Ok"), SEND_REPLY_TIMEOUT);
}

template<class T>
void udpSendTest(TcpConnection& conn, const T& ref)
{
	T data(ref.size(), 0);
	UdpTransceiver trans(SocketAddress(IpAddress(127, 0, 0, 1), 10002));
	SocketAddress source;
	trans.doRcve(source, buffer(data), RCVE_DATA_TIMEOUT);
	checkResult(conn, ref, data);
}

void udp_send1(TcpConnection& conn)
{
	NetTestData ref;
	udpSendTest(conn, ref.dataset1);
}

void udp_send2(TcpConnection& conn)
{
	NetTestData ref;
	udpSendTest(conn, ref.dataset1);
}

void udp_recv1(TcpConnection& conn)
{
	NetTestData data;
	UdpTransceiver trans;
	trans.doSend(SocketAddress(IpAddress(127, 0, 0, 1), 10003), buffer(data.dataset3));
	conn.doSendValue(std::string("Ok"), SEND_REPLY_TIMEOUT);
}

void udp_recv_timeout(TcpConnection& conn)
{
	mdaq::mt::SysClock::sleep_us(500000);
	conn.doSendValue(std::string("Ok"), SEND_REPLY_TIMEOUT);
}

#define REGISTER_TEST(TESTNAME) testMap.insert(std::make_pair(std::string( #TESTNAME ), TESTNAME))


int main()
{
	typedef std::map<std::string, void (*)(TcpConnection&)> TestMap;
	TestMap testMap;
	REGISTER_TEST(tcp_send1);
	REGISTER_TEST(tcp_send2);
	REGISTER_TEST(tcp_send4);
	REGISTER_TEST(tcp_send5);
	REGISTER_TEST(tcp_send6);
	REGISTER_TEST(tcp_send_timeout);
	REGISTER_TEST(udp_send1);
	REGISTER_TEST(udp_send2);
	REGISTER_TEST(udp_recv1);
	REGISTER_TEST(udp_recv_timeout);

	try
	{
		bool keepAlive = true;
		do
		{
			TcpAcceptor acceptor(SocketAddress(IpAddress(127, 0, 0, 1), 10001));
			std::auto_ptr<TcpConnection> conn = acceptor.acceptConnection(ACCEPT_TIMEOUT);
			std::string testName;
			conn->doRcveValue(testName, RCVE_DATA_TIMEOUT);
			if(testName != "quit")
			{
				TestMap::const_iterator it = testMap.find(testName);
				if(it != testMap.end())
				{
					(*it->second)(*conn);
				}
				else
				{
					LOG_WARN() << "Unknown test " << testName;
					conn->doSendValue(std::string("Unknown test"), SEND_REPLY_TIMEOUT);
				}
			}
			else
			{
				keepAlive = false;
			}
		}
		while(keepAlive);
	}
	catch(const std::exception& e)
	{
		LOG_ERROR() << "Exception in server: " << e.what();
	}
	catch(...)
	{
		LOG_ERROR() << "Exception in server: " << "unknown exception";
	}
	return 0;
}


