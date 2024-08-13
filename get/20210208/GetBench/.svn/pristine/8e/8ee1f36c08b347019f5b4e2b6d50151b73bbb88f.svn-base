/**
 * @file FakeFrameSender.cpp
 * @date 19 nov. 2012
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

#include "get/daq/FakeFrameSender.h"
#include "utl/Logging.h"
#include "mfm/Exception.h"
#include "mfm/Frame.h"
#include <iomanip>
#include <iostream>

namespace get
{
namespace daq
{
//_________________________________________________________________________________________________
FakeFrameSender::FakeFrameSender()
	: bytesSent(0), framesSent(0)
{
}
//_________________________________________________________________________________________________
FakeFrameSender::~FakeFrameSender()
{
}
//_________________________________________________________________________________________________
void FakeFrameSender::sendData(const uint8_t* buffer, uint32_t nBytes)
{
	const mfm::Byte* chunkBeg =  reinterpret_cast<const mfm::Byte*>(buffer);
	const mfm::Byte* chunkEnd =  reinterpret_cast<const mfm::Byte*>(buffer + nBytes);

	try
	{
		FrameBuilder::addDataChunk(chunkBeg, chunkEnd);
		bytesSent += nBytes;
		return;
	}
	catch (const mfm::Exception& e)
	{
		// This is probably due to a data format error: dump some of it for debugging:
		LOG_ERROR() << e.what();
		LOG_ERROR() << "Dumping data bytes in hexadecimal (but no more than 8 lines, i.e. 256 Bytes): ";
		std::ios::fmtflags flags(std::cerr.flags()); // saving format flags
		std::cerr << std::hex << std::right << std::setfill('0');
		for (size_t i = 0; (i < nBytes) and (i < 8*32); ++i)
		{
			if (i % 32 == 0) std::cerr << '\n' << std::setw(3) << i << ": ";
			if (i % 8 == 0) std::cerr << ' ';
			std::cerr << std::setw(2) << (int) buffer[i] << ' ' ;
		}
		std::cerr << std::endl;
		std::cerr.flags(flags); // restoring format flags

	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
	}
	// Reset frame builder
	LOG_WARN() << "Resetting frame builder...";
	FrameBuilder::reset();
}
//_________________________________________________________________________________________________
void FakeFrameSender::processFrame(mfm::Frame & frame)
{
	framesSent++;
}
//_________________________________________________________________________________________________
void FakeFrameSender::getStatistics(uint64_t & byteCount, uint64_t & frameCount) const
{
	byteCount =  bytesSent;
	frameCount = framesSent;
}
//_________________________________________________________________________________________________
} /* namespace daq */
} /* namespace get */
