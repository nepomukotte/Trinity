/* =====================================================================================================================
 * CmdLineArgs.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::utl::CmdLineArgs
 * Created on: Jan 11, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the MDaq data acquisition framework @ CEA-Irfu.
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
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#include "CmdLineArgs.h"
#include <boost/tokenizer.hpp>
#include <string.h>

namespace mdaq
{
namespace utl
{

CmdLineArgs CmdLineArgs::build(int argc, char** argv)
{
	CmdLineArgs args;
	for (int i=0; i < argc; ++i)
	{
		args.push_back(argv[i]);
	}
	return args;
}

/**
 * Constructs a CmdLineArgs object by copying a vector of argument strings.
 * @param args Arguments.
 */
CmdLineArgs::CmdLineArgs(const std::vector< std::string > & args)
	: std::vector< std::string >(args), argv_(0)
{
	;
}

void CmdLineArgs::destroy_argv()
{
	if (argv_ != 0)
	{
		for (char* c_str = argv_[0]; *c_str != 0; ++c_str)
		{
			delete [] c_str;
		}
	}
	delete [] argv_;
	argv_ = 0;
}

CmdLineArgs CmdLineArgs::build(const char* funcName, const char* funcArgs)
{
	std::string funcArgsStr(funcArgs);
	boost::char_separator<char> sep(" "); // only whitespace separator
	boost::tokenizer< boost::char_separator<char> > tokens(funcArgsStr, sep); // tokenizer instantiated
	CmdLineArgs args;
	args.push_back(funcName);
	for(boost::tokenizer< boost::char_separator<char> >::iterator arg = tokens.begin(); arg != tokens.end(); ++arg)
	{
		args.addArg(*arg);
	}
	return args;
}

CmdLineArgs& CmdLineArgs::erase(size_type i)
{
	if (i < size())
	{
		iterator it = begin();
		it += i;
		std::vector<std::string>::erase(it);
	}
	return *this;
}

char** CmdLineArgs::argv()
{
	destroy_argv();

	argv_ = new char* [1+size()];
	argv_[size()] = 0;
	for (size_type i=0; i < size(); ++i)
	{
		argv_[i] = new char [1+at(i).size()];
		argv_[i] = strcpy(argv_[i], at(i).c_str());
	}
	return argv_;
}

} // namespace utl
} // namespace mdaq
