/* =====================================================================================================================
 * Cmd.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::cli::Cmd
 * class mdaq::cli::Cmd::Engine
 * Created on: Jan 5, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
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

#ifndef mdaq_cli_Cmd_INCLUDED
#define mdaq_cli_Cmd_INCLUDED

#include "mdaq/utl/net/SocketAddress.h"
#include "UserQuitting.h"
#include <Ice/Ice.h>
#include <Ice/Application.h>
#include <boost/algorithm/string/case_conv.hpp>
#include <map>
#include <memory>
#include <fstream>
#include <readline/readline.h>

namespace mdaq
{
namespace cli
{
class Cmd
{
public:
	/// This class extends a map by adding cleanup of the contained heap-allocated pointers.
	struct Map : public std::map<std::string, Cmd*>
	{
		virtual ~Map() { for (iterator it=begin(); it != end(); ++it) delete it->second; }
	};
	typedef std::map<std::string, std::string> AliasDict;

	/// This class implements the command engine.
	class Engine : virtual public Ice::Application
	{
	public:
		Engine(const std::string& name);
		virtual ~Engine();

		virtual int run(int argc, char*argv[]); // Inherited from Ice::Application: to be implemented by derived class
		void execNextCmd(); // Body of the run loop

		Cmd* searchCmd(const std::string& cmdName);
		Cmd& findCmd(const std::string& cmdName);
		Cmd& nextCmd(std::istream& args);

		template <class CmdClass>
		void subscribe();

		void addAlias(const std::string & cmdAlias, const std::string & cmdName);
		Cmd* searchAlias(const std::string& cmdAlias);

		void help(std::string& cmdName);
		void help();

	public:
		const std::string name;
		virtual std::string prompt() { return idStr() + " > "; }
		virtual std::string idStr() { return name; }

	protected:
		Cmd::Map cmdMap;
		Cmd::AliasDict cmdAliases;
		::mdaq::utl::net::SocketAddress ip;

	private:
		void initHistPath();
		std::string histPath;
		static const int maxHist = 200;
	};

	/**
	 * Exception safe command line
	 */
	class Line
	{
	public:
		Line(const std::string& prompt) { buffer = readline(prompt.c_str()); }
		~Line() { free(buffer); }
		const std::string toString() { return buffer == 0 ? "" : buffer; }
		operator std::string () { return toString(); }

	private:
		char* buffer;
	};

public:
	virtual ~Cmd() {};
	void help(void);
	void completeHelp(void);
	std::string usage() { return name + ' ' + argsDescription; }
	virtual void exec(std::istream& cmdArgs) = 0;
	static std::string to_lower(const std::string& str) { std::string s(str); boost::algorithm::to_lower(s); return s; }

protected:
	Cmd(Engine& engine, const std::string& name, const std::string& argsDescription, const std::string& comment)
		: engine(engine), name(to_lower(name)), argsDescription(argsDescription), comment(comment) {}

protected:
	Engine& engine;

public:
	const std::string name;
	const std::string argsDescription;
	const std::string comment;
};

template<class CmdClass>
void Cmd::Engine::subscribe()
{
	Cmd* cmd = new CmdClass(dynamic_cast<typename CmdClass::Engine&>(*this));
	Cmd::Map::iterator it = cmdMap.find(cmd->name);
	if (it != cmdMap.end())
	{
		delete it->second;
	}
	cmdMap[cmd->name] = cmd;
}

} //cli
} //mdaq

#endif /* mdaq_cli_Cmd_INCLUDED */
