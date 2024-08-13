/**
 * @file ByteStorage.cpp
 * @date 31 oct. 2012
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

#include "ByteStorage.h"
#include "get/utl/Time.h"
#include "utl/Logging.h"
#include <iomanip>

namespace get
{
namespace daq
{

ByteStorage::ByteStorage() : fileName(), maxFileSize_MiB(1024u)
{
	LOG_DEBUG() << "Creating ByteStorage " << std::hex << this << std::dec;
	createNewFile("ByteStorage");
}

ByteStorage::~ByteStorage()
{
	closeFile();
	LOG_DEBUG() << "Destroying ByteStorage";
}

void ByteStorage::resetData()
{
}

/**
 * Adds data to frame which is currently being reconstructed.
 */
void ByteStorage::processData(const uint8_t* buffer, uint32_t nBytes)
{
	try
	{
		LOG_DEBUG() << "Storing chunk (" << nBytes << " B)";
		openFile();
		const char* data =  reinterpret_cast<const char*>(buffer);
		dataFile.write(data, nBytes);
		// Close file after recording its size
		const size_t fileSize_B = dataFile.tellp();
		closeFile();

		// Create new file is maximum file size has been reached
		if (fileSize_B > maxFileSize_MiB*0x100000)
			createNewFile("ByteStorage");
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
	}
}

void ByteStorage::openFile()
{
	if (not dataFile.is_open())
	{
		dataFile.open(fileName.c_str(), std::ios::out | std::ios::app);
		if (dataFile.is_open())
		{
			LOG_DEBUG() << "Appending data to file: " << fileName;
		}
		else
		{
			LOG_ERROR() << "Could not open file : " << fileName;
		}
	}
}

void ByteStorage::closeFile()
{
	if (dataFile.is_open())
		dataFile.close();
}

void ByteStorage::createNewFile(const std::string & prefix)
{
	closeFile();
	std::ostringstream oss;
	oss << prefix << '_' << get::utl::buildTimeStamp() << ".graw";
	fileName = oss.str();
	dataFile.open(fileName.c_str(), std::ios::out | std::ios::trunc);
	closeFile();
	LOG_INFO() << "New file: " << fileName;
}

} /* namespace daq */
} /* namespace get */
