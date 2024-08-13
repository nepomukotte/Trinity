/* =====================================================================================================================
 * DirectStorage.h
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

#ifndef get_daq_FrameStorage_h_INCLUDED
#define get_daq_FrameStorage_h_INCLUDED

#include "mfm/FrameBuilder.h"
#include "mdaq/daq/DataProcessorCore.h"
#include <cstddef>
#include <map>
#include <memory>
#include <fstream>
#include <string>
#include <utility>

namespace get
{
namespace daq
{

class FrameStorage : public mdaq::daq::DataProcessorCore, public virtual mfm::FrameBuilder
{
public:
	FrameStorage();
	virtual ~FrameStorage();
	virtual void resetData(); ///< Overrides method inherited from mdaq::daq::DataProcessorCore
	virtual void processData(const uint8_t* buffer, uint32_t nBytes); ///< Overrides method inherited from mdaq::daq::DataProcessorCore
	virtual void processFrame(mfm::Frame& frame); ///< Implements method inherited from mfm::FrameBuilder
	virtual void processHeader(const mfm::PrimaryHeader & header);  ///< Overrides method inherited from mfm::FrameBuilder
	typedef std::pair< size_t, size_t > SourceUID; ///< CoBo/MuTanT index, AsAd index
	std::ofstream & getFile(const SourceUID & sourceIdx, const size_t & eventIdx);
	std::ofstream & openFile(const SourceUID & sourceIdx);
	void setOutputDirectory(const std::string & dir);
	void closeFile(const SourceUID & sourceIdx);
	void createNewFile(const std::string & prefix, const SourceUID & asadIdx, bool newTimeStamp=false);
	/**
	 * Structure grouping information pertaining to storage of frames for a specific AsAd board or MuTanT module.
	 */
	struct StorageInfo
	{
		uint32_t previousEventIdx; ///< Event ID of last event.
		std::ofstream dataFile; ///< File for storage.
		std::string fileName; ///< File name.
		size_t fileIndex; ///< Index of file.
		StorageInfo();
		StorageInfo(const StorageInfo &);
	};
private:
	FrameStorage(const FrameStorage& r);
	FrameStorage& operator=(const FrameStorage& r);
private:
	std::map < SourceUID, StorageInfo > storages;
	std::string filePrefix;
	std::string fileTimestamp;
	size_t maxFileSize_MiB;
	std::string fileDirectory; ///< Path of directory to start the frames to.
};

} // namespace daq
} // namespace get

#endif
