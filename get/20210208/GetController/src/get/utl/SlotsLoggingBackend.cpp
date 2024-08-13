/**
 * @file SlotsLoggingBackend.cpp
 * @date 18 oct. 2012
 * @author fchateau
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Frederic Chateau, Patrick Sizun
 *
 * This file is part of the GetController software project.
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

#include "get/utl/SlotsLoggingBackend.h"
#include <cassert>
#include <iomanip>

namespace get
{
namespace utl
{

SlotsLoggingBackend::SlotsLoggingBackend() : errorMode(0), mutex(QMutex::Recursive)
{
}

SlotsLoggingBackend::~SlotsLoggingBackend()
{
	QObject::disconnect(this, 0, 0, 0);
}

std::ostream& SlotsLoggingBackend::beginLog(::utl::SeverityLevel level, const char* /*source*/, unsigned /*line*/)
{
	static const char* prompt[::utl::fatal+1] = {"==", "==", "!!", "**", "##"};
	static const char* color[::utl::fatal+1] = {"black", "blue", "orange", "red", "red"};
	assert(level <= ::utl::fatal);
	mutex.lock();
	errorMode = (level >= ::utl::error);
	stream.str(std::string());
	stream << std::resetiosflags(std::ios_base::showbase) << std::setbase(0);
	stream << "<font color=" << color[level] << '>' << prompt[level] << ' ';
	return stream;
}

void SlotsLoggingBackend::endLog()
{
	stream << "</font>";
	if (errorMode)
		emit error(QString::fromStdString(stream.str()));
	else
		emit log(QString::fromStdString(stream.str()));
	mutex.unlock();
}

} /* namespace utl */
} /* namespace get */
