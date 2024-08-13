/* =============================================================================
 * DataBloc.h
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
#ifndef get_DataBloc_h_INCLUDED
#define get_DataBloc_h_INCLUDED

#include <QObject>
#include <QMutex>
#include <QMutexLocker>
#include <QWaitCondition>
#include <boost/cstdint.hpp>
#include <vector>
#include <deque>
#include <map>

/**
 * Thread-safe frame raw data circular buffer.
 */
class DataBloc : public QObject
{
	 Q_OBJECT /*!< Belong to Qt signal/slot interface*/
public:
	 typedef std::vector< char > RawFrame;
	 typedef std::deque<RawFrame*> RawQueue;
	 typedef std::multimap<size_t, RawFrame*> RawPool;
	 static const size_t MEM_BUDGET_MiB = 64u; // Maximum allocated data memory in megabytes.
public:
	DataBloc();
	virtual ~DataBloc();
	void clear();
public:
	RawFrame* prepareRawFrame(size_t capacity_B);
	void pushFrame(RawFrame* rawFrame);
	RawFrame* waitNextFrame(unsigned long timeout_ms = ULONG_MAX);
	void popFrame(RawFrame* rawFrame);
	size_t byteCount() const {  return byteCount_; }
	size_t frameCount() const { QMutexLocker locker(&mutex); return rawQueue.size(); }
	size_t allocCount() const {  return allocCount_; }
	size_t droppedFrames() const {  return droppedFrames_; }
	size_t droppedBytes() const {  return droppedBytes_; }
private:
	void drop(size_t nBytes);
private:
	RawQueue rawQueue;
	size_t byteCount_;
	size_t allocCount_;
	size_t droppedFrames_;
	size_t droppedBytes_;
	RawPool rawPool;
	mutable QMutex mutex;
	QWaitCondition notEmpty;
};

#endif /* get_DataBloc_h_INCLUDED */
