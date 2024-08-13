/* =============================================================================
 * BuildInfo.cpp
 * -----------------------------------------------------------------------------
 * class utl::BuildInfo
 * Created on: 7 août 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This software is part of the Utilities project.
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

#include "config.h"
#include "genbuildinfo.h"
#include "BuildInfo.h"
#include <iostream>

namespace utl
{

BuildInfo* BuildInfo::first_ = 0;

/**
 * Constructor.
 */
BuildInfo::BuildInfo(const char* module, const char* date, const char* time, const char* compiler, const char* revision)
	: MODULE(module),
	  DATE(date),
	  TIME(time),
	  COMPILER_VERSION(compiler),
	  REVISION(revision),
	  next_(0),
	  prev_(0)
{
	if(first_ == 0)
	{
		first_ = this;
		next_ = prev_ = this;
	}
	else
	{
		BuildInfo* prev = first_->prev_;
		BuildInfo* next = first_;
		this->next_ = next;
		this->prev_ = prev;
		next->prev_ = this;
		prev->next_ = this;
	}
}

BuildInfo::~BuildInfo()
{
	if(next_ == prev_)
	{
		first_ = 0;
	}
	else
	{
		this->next_->prev_ = this->prev_;
		this->prev_->next_ = this->next_;
	}
	this->next_ = 0;
	this->prev_ = 0;
}

void BuildInfo::printAll(std::ostream& s)
{
	BuildInfo* endInfo = BuildInfo::first()->prev();
	BuildInfo* info = endInfo;
	do
	{
		if(!info) break;
		s << *info << std::endl;
		info = info->prev();
	}
	while(info != endInfo);
}

std::ostream& operator<<(std::ostream& s, const BuildInfo& b)
{
	s << "Module \'" << b.MODULE << "\' build time: " << b.DATE << ' ' << b.TIME
	  << ", compiler: \'" << b.COMPILER_VERSION << "\', sources revision: " << b.REVISION;
	return s;
}

} /* namespace utl */


// including build informations for the Utilities library.
BUILD_INFO

