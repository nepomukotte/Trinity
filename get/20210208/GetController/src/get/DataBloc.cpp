/* =============================================================================
 * DataBloc.cpp
 * -----------------------------------------------------------------------------
 * class get::DataBloc
 * Created on: 29 nov. 2010 at Irfu/Sedi/Ldef, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Druillole (frederic.druillole@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the GetController project dedicated to the GET
 * data acquisition test bench. "GET" stands for General Electronics for TPC.
 * http://www-actar-get.cea.fr
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

#include "get/DataBloc.h"
#include "utl/Logging.h"
#include <iomanip>
#include <boost/cstdint.hpp>

#define DB_LOG(expr)
#define ERR_LOG(expr)

DataBloc::DataBloc() : mutex(QMutex::NonRecursive)
{
	LOG_DEBUG() << "Creating queue for frames";
	clear();
}

void DataBloc::clear()
{
	LOG_WARN() << "Resetting data buffer.";
	QMutexLocker locker(&mutex);

	for (RawPool::iterator it = rawPool.begin(); it != rawPool.end(); ++it)
	{
		delete it->second;
	}
	rawPool.clear();

	while (rawQueue.size() > 0)
	{
		delete rawQueue.front();
		rawQueue.pop_front();
	}

	byteCount_ = 0u;
	allocCount_ = 0u;
	droppedFrames_ = 0u;
	droppedBytes_ = 0u;
}

DataBloc::~DataBloc()
{
	clear();
}

DataBloc::RawFrame* DataBloc::prepareRawFrame(size_t capacity_B)
{
	QMutexLocker locker(&mutex);
	RawPool::iterator it_rawFrame = rawPool.lower_bound(capacity_B);
	RawFrame* rawFrame = 0;
	if (it_rawFrame == rawPool.end())
	{
		if (allocCount_ + capacity_B <= MEM_BUDGET_MiB*0x100000u)
		{
			rawFrame = new RawFrame;
			rawFrame->reserve(capacity_B);
			allocCount_ += capacity_B;

			LOG_DEBUG() << "Allocating " << std::dec << capacity_B << " B array to store frame ("
					<< std::setfill('0') << std::setw(5) << std::setprecision(5) << (allocCount_*1.0/0x100000) << " MiB allocated)";
		}
		else
		{
			drop(capacity_B); // Data loss
			rawFrame = 0;
		}
	}
	else
	{
		rawFrame = it_rawFrame->second;
		LOG_DEBUG() << "Reusing array of " << rawFrame->capacity() << " B.";
		rawPool.erase(it_rawFrame);
		rawFrame->clear();
	}
	return rawFrame;
}

void DataBloc::pushFrame(RawFrame* rawFrame)
{
	if (rawFrame != 0)
	{
		QMutexLocker locker(&mutex);
		rawQueue.push_back(rawFrame);
		byteCount_ += rawFrame->size();
		notEmpty.wakeAll();
	}
}

DataBloc::RawFrame* DataBloc::waitNextFrame(unsigned long timeout_ms)
{
	QMutexLocker locker(&mutex);
	while (rawQueue.size() == 0u)
	{
		if (not notEmpty.wait(&mutex, timeout_ms))
		{
			return 0;
		}
	}
	RawFrame* rawFrame = rawQueue.front();
	rawQueue.pop_front();
	byteCount_ -= rawFrame->size();
	return rawFrame;
}

void DataBloc::popFrame(RawFrame* rawFrame)
{
	if (rawFrame != 0)
	{
		rawFrame->clear();
		QMutexLocker locker(&mutex);
		rawPool.insert(std::make_pair(rawFrame->capacity(), rawFrame));
	}
}

void DataBloc::drop(size_t nBytes)
{
	LOG_DEBUG() << "Dropping frame of " << std::dec << nBytes << " Bytes";
	++droppedFrames_;
	droppedBytes_ += nBytes;
}
