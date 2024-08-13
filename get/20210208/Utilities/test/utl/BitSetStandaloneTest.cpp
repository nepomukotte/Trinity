/*
 * BitSetSmallTest.cpp
 *
 *  Created on: Apr 11, 2016
 *      Author: anvar
 */
#include <utl/BitSet.hpp>
#include <vector>
#include <iostream>
#include <iomanip>
#include <boost/cstdint.hpp>

using namespace std;

int main()
{
	vector<uint16_t> buf(16, 0xface);
	for (size_t i=0; i < buf.size(); ++i) buf[i] = (i+1) | ((i+2) << 8);
	cout << hex << setfill('0');
	for (unsigned i=0; i < buf.size(); ++i) cout << setw(4) << buf[i] << ", "; cout << endl;

	utl::BitSet<123, uint16_t> bs(&buf[0]);
	cout << bs << endl;

	cout << setw(2) << bs.getField<int>(0, 8) << endl;
	bs.setField(10, 8, 0x11);
	cout << bs << endl;

	for (unsigned i=0; i < buf.size(); ++i) cout << setw(4) << buf[i] << ", "; cout << endl;
	cout << hex << setw(2) << bs.getField<int>(0, 8) << endl;
	return 0;
}

