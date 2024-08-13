/* =====================================================================================================================
 * DirectStorage.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class get::daq::DirectStorage
 * Created on: Apr 4, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * (c) Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 *               Patrick Sizun (patrick.sizun@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GetBench package implementing the GET-specific software modules for the MDaq framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
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
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#include "FrameStorage.h"
#include "get/utl/Time.h"
#include "mfm/Serializer.h"
#include "mfm/Frame.h"
#include "mfm/FrameDictionary.h"
#include "utl/Logging.h"
#include <iomanip>
#include <sstream>
#include <utility>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

namespace get
{
namespace daq
{

FrameStorage::FrameStorage() : maxFileSize_MiB(1024u), fileDirectory(".")
{
	LOG_DEBUG() << "Creating FrameStorage " << std::hex << this << std::dec;
	// Load Cobo format
	//mfm::FrameDictionary::instance().addFormats("CoboFormats.xcfg");
}

FrameStorage::StorageInfo::StorageInfo() : previousEventIdx(1), fileName(), fileIndex(0)
{
}

/**
 *  Copy constructor.
 */
FrameStorage::StorageInfo::StorageInfo(const StorageInfo & r)
 :
	previousEventIdx(r.previousEventIdx),
	fileName(r.fileName),
	fileIndex(r.fileIndex)
{
}

FrameStorage::~FrameStorage()
{
	typedef std::map<SourceUID, FrameStorage::StorageInfo>::iterator StorageInfoIter;
	for (StorageInfoIter it = storages.begin(); it != storages.end(); ++it)
	{
		closeFile(it->first);
	}
	LOG_DEBUG() << "Destroying FrameStorage";
}

/**
 * Sets path of directory to save frames to.
 */
void FrameStorage::setOutputDirectory(const std::string & dir)
{
	fs::path outputPath(dir);
	if (not fs::exists(outputPath))
	{
		LOG_ERROR() << "Directory '" << dir << "' does not exist!";
		return;
	}
	fileDirectory = dir;
}

void FrameStorage::resetData()
{
	reset(); // Discard pending data frame chunks
}

/**
 * Adds data to frame which is currently being reconstructed.
 */
void FrameStorage::processData(const uint8_t* buffer, uint32_t nBytes)
{
	const mfm::Byte* chunkBeg =  reinterpret_cast<const mfm::Byte*>(buffer);
	const mfm::Byte* chunkEnd =  reinterpret_cast<const mfm::Byte*>(buffer + nBytes);

	try
	{
		FrameBuilder::addDataChunk(chunkBeg, chunkEnd);
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

/**
 * Writes reconstructed frame to disk, creating a new file if necessary.
 */
void FrameStorage::processFrame(mfm::Frame& frame)
{
	// Get file
	uint32_t coboIdx, asadIdx, eventIdx;
	SourceUID sourceUID;

	try
	{
		coboIdx = frame.headerField("coboIdx").value< uint8_t >();
		asadIdx = frame.headerField("asadIdx").value< uint8_t >();
		sourceUID = std::make_pair<>(coboIdx, asadIdx);
		eventIdx = frame.headerField("eventIdx").value<uint32_t>();
		LOG_DEBUG() << "Storing AsAd no. " << asadIdx << " frame for event " << std::dec << eventIdx << " (" << frame.header().frameSize_B() << " B)";
	}
	catch (const mfm::Exception &)
	{
		// Check if frame is of MuTanT type 0
		if (frame.header().frameType() == 0x8)
		{
			sourceUID = std::make_pair<>(100, 0);
			eventIdx = frame.headerField(14u, 4u).value< uint32_t >();
			LOG_DEBUG() << "Storing MuTanT frame for event " << std::dec << eventIdx << " (" << frame.header().frameSize_B() << " B)";
		}
		else
		{
			LOG_WARN() << "Dumping frame of type " << std::hex << std::showbase << frame.header().frameType() << std::dec;
			return;
		}
	}

	std::ofstream & file = getFile(sourceUID, eventIdx);

	// Store frame in file

	frame.write(file);

	// Deal with write errors (no space left on device)
	if (not file.good())
	{
		LOG_ERROR() << "Error writing frame to file!";
		closeFile(sourceUID);
		return;
	}

	// Close file after recording its size
	const size_t fileSize_B = file.tellp();
	closeFile(sourceUID);

	// Create new file if maximum file size has been reached
	if (fileSize_B > maxFileSize_MiB*0x100000)
	{
		createNewFile("", sourceUID, false);
	}
}

void FrameStorage::processHeader(const mfm::PrimaryHeader & header)
{
	LOG_DEBUG() << "Building frame of " << header.frameSize_B() << " B";
}

/**
 * Returns open file for storing frame, creating a new file if necessary.
 */
std::ofstream & FrameStorage::getFile(const SourceUID & sourceIdx, const size_t & eventIdx)
{
	StorageInfo & storage = storages[sourceIdx];

	// Create new file if eventIdx has been reset
	if (eventIdx <= storage.previousEventIdx or storage.fileName == "")
	//if (fileName == "")
	{
		createNewFile("", sourceIdx, true);
	}
	storage.previousEventIdx = eventIdx;

	return openFile(sourceIdx);
}

/**
 * Opens existing file for storing frames for an AsAd board.
 * @param sourceIdx ID of AsAd board or MuTanT module.
 */
std::ofstream & FrameStorage::openFile(const SourceUID & sourceIdx)
{
	StorageInfo & storage = storages[sourceIdx];

	if (not storage.dataFile.is_open())
	{
		storage.dataFile.clear();
		storage.dataFile.open(storage.fileName.c_str(), std::ios::binary | std::ios::out | std::ios::app);
		if (storage.dataFile.is_open())
		{
			LOG_DEBUG() << "Appending data to file: " << storage.fileName;
		}
		else
		{
			LOG_ERROR() << "Could not open file : " << storage.fileName;
		}
	}
	return storage.dataFile;
}

/**
 * Closes file for storing frames for an AsAd board.
 * @param asadIdx ID of AsAd board or MuTanT module.
 */
void FrameStorage::closeFile(const SourceUID & sourceIdx)
{
	StorageInfo & storage = storages[sourceIdx];
	if (storage.dataFile.is_open())
		storage.dataFile.close();
}

/**
 * Creates a new file for storages of frames of a specific AsAd board.
 * @param prefix Filename prefix.
 * @param sourceIdx ID of AsAd board or MuTanT module.
 * @param newTimeStamp Whether to create a new timestamp for the filename.
 */
void FrameStorage::createNewFile(const std::string & prefix, const SourceUID & sourceIdx, bool newTimeStamp)
{
	closeFile(sourceIdx);

	StorageInfo & storage = storages[sourceIdx];

	// File timestamp and index
	if (newTimeStamp or fileTimestamp.empty())
	{
		fileTimestamp = get::utl::buildTimeStamp();
		storage.fileIndex = 0;
	}
	else
	{
		storage.fileIndex++;
	}
	// File prefix
	if (not prefix.empty())
	{
		filePrefix = prefix + '_';
	}
	// Filename
	std::ostringstream oss;
	oss << fileDirectory << '/' << filePrefix;
	if (sourceIdx.first < 100)
	{
		oss << "CoBo" << sourceIdx.first << "_AsAd" << sourceIdx.second;
	}
	else
	{
		oss << "MuTanT";
	}
	oss << '_' << fileTimestamp << '_' << std::setfill('0') << std::setw(4) << storage.fileIndex << ".graw";
	storage.fileName = oss.str();

	storage.dataFile.clear();
	storage.dataFile.open(storage.fileName.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
	// Deal with disk with no space left
	if (storage.dataFile.fail())
	{
		LOG_ERROR() << "Could not create file " << storage.fileName;
	}
	closeFile(sourceIdx);
	LOG_INFO() << "New file: " << storage.fileName;
}

} // namespace daq
} // namespace get

