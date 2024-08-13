/* =====================================================================================================================
 * Cmd.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq:cli::Cmd
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

#include "Cmd.h"
#include "mdaq/utl/CmdLineArgs.h"
#include "mdaq/utl/CmdException.h"
#include "UserQuitting.h"
#include <readline/history.h>
#include <sstream>
#include <iomanip>
#include <string>

namespace mdaq
{
namespace cli
{

/**
 * This command class allows the execution of command files.
 */
struct Exec : public Cmd
{
	Exec(Engine& engine) : Cmd(engine, "exec", "<filename>", "Executes a command file") {}
	void exec(std::istream& cmdArgs)
	{
		std::string execFileName;
		cmdArgs >> execFileName;

		// Opens the file and start reading and executing commands line by line
		std::ifstream cmdFile(execFileName.c_str());
		if (! cmdFile)
		{
			throw utl::CmdException(std::string("could not open command file '") + execFileName + '\'');
		}
		std::string cmdLine;
		while (std::getline(cmdFile, cmdLine))
		{
			if (not cmdLine.empty())
			{
				const std::string::size_type pos = cmdLine.find_first_not_of(" \t");
				if (std::string::npos != pos)
					cmdLine = cmdLine.substr(pos); // Skip possible leading spaces
				else
					cmdLine.clear();
			}
			if (not cmdLine.empty())
			{
				cmdLine = cmdLine.substr(0, cmdLine.find_first_of("#")); // Ignore everything after # character
			}
			if (cmdLine != "")
			{
				std::stringstream cmdLineS(cmdLine);
				std::cout << engine.prompt() << cmdLine << std::endl;
				engine.nextCmd(cmdLineS).exec(cmdLineS);
			}
		}
		cmdFile.close();
	}
};

/**
 * This command is the NOP command: it does nothing.
 */
struct Nop : public Cmd
{
	Nop(Engine& engine) : Cmd(engine, "nop", "", "Executes nothing") {}
	void exec(std::istream&) {}
};

/**
 * This command executes a command the specified number of times.
 */
struct Loop : public Cmd
{
	Loop(Engine& engine)
		: Cmd(engine,
			"loop",
			"<count> <command> [arguments]",
			"Loops count times over the execution of specified command every half-second"
			) {}

	void exec(std::istream& cmdArgs)
	{
		int count;
		cmdArgs >> count;
		if (cmdArgs.fail())
			throw usage();

		std::string cmdLine;
		std::getline(cmdArgs, cmdLine);
		while (count-- > 0)
		{
			std::istringstream is(cmdLine);
			engine.nextCmd(is).exec(is);
			IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(500));
		}
	}
};

/**
 * This command indefinitely executes a command until a key is pressed.
 */
struct LoopKey : public Cmd, public IceUtil::Thread
{
	volatile bool stopLoop;
	virtual void run() {
		stopLoop = false;
		std::cin.get();
		stopLoop = true;
	};

	LoopKey(Engine& engine)
		: Cmd(engine,
			"loopkey",
			"<command> [arguments]",
			"Loops the execution of specified command every half-second until any <RETURN> key is pressed"
			) {}

	void exec(std::istream& cmdArgs)
	{
		std::string cmdLine;
		std::getline(cmdArgs, cmdLine);

		std::cout << "To stop the loop, press <RETURN>" << std::endl;
		stopLoop = false;
		start();
		while (not stopLoop)
		{
			std::istringstream is(cmdLine);
			engine.nextCmd(is).exec(is);
			IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(500));
		}
	}
};

/**
 * This command is useful for scripting as it allows script execution to sleep for specified amount of time.
 */
struct Sleep : public Cmd
{
	Sleep(Engine& engine) : Cmd(engine, "sleep", "<seconds>", "Sleeps specified number of seconds") {}
	void exec(std::istream& cmdArgs)
	{
		uint64_t nSeconds;
		cmdArgs >> nSeconds;
		if (cmdArgs.fail())
			throw usage();

		IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(nSeconds*1000));
	}
};

/**
 * This command prints a command help.
 */
struct Help : public Cmd
{
	Help(Engine& engine) : Cmd(engine, "help", "[command]", "Help on a specific command") {}
	void exec(std::istream& args)
	{
		std::string cmdName;
		args >> cmdName;
		engine.help(cmdName);
	}
};

/**
 * This command quits the command line interpreter.
 */
struct Quit : public Cmd
{
	Quit(Engine& engine) : Cmd(engine, "q", "", "Exits the program") {}
	void exec(std::istream&)
	{
		throw UserQuitting();
	}
};

Cmd::Engine::Engine(const std::string& name) : name(name)
{
	initHistPath();

	// enable command line history:
	using_history();
	stifle_history(maxHist);

	// Read out previous history
	read_history(histPath.c_str());

	// Default commands
	subscribe<Nop>();
	subscribe<Exec>();
	subscribe<Loop>();
	subscribe<LoopKey>();
	subscribe<Sleep>();
	subscribe<Help>();
	subscribe<Quit>();
}

Cmd::Engine::~Engine()
{
	if (write_history(histPath.c_str()))
		perror("Saving history: ");
}

int Cmd::Engine::run(int argc, char*argv[])
{
	utl::CmdLineArgs args = utl::CmdLineArgs::build(argc, argv);
	if (args.size() > 1)
	{
		ip.fromString(args[1]);
	}

	help();
	while (true) // Interpretation of the line commands
	{
		try
		{
			execNextCmd();
		}
		catch (UserQuitting& e)
		{
			return 0; // Exit command loop.
		}
		catch (const utl::CmdException & e)
		{
			std::cerr << "*** Error: " << e.reason << std::endl;
		}
		catch (const Ice::Exception & e)
		{
			std::cerr << "*** " << e.what() << std::endl;
		}
		catch (const std::exception& e)
		{
			std::cerr << "*** " << e.what() << std::endl;
		}
		catch (const std::string& msg)
		{
			std::cerr << "*** " << msg << std::endl;
		}
		catch (const char* msg)
		{
			std::cerr << "*** " << msg << std::endl;
		}
		catch(...)
		{
			std::cerr << "*** Unknown error" << std::endl;
		}
	}
	return 1;
}

void Cmd::Engine::execNextCmd()
{
	std::string cmdLine = Cmd::Line(prompt());
	if (cmdLine != "")
	{
		add_history(cmdLine.c_str()); //save line in history
		std::istringstream cmdLineS(cmdLine);
		nextCmd(cmdLineS).exec(cmdLineS);
	}
}

Cmd* Cmd::Engine::searchCmd(const std::string& cmdName)
{
	Cmd::Map::iterator it = cmdMap.find(to_lower(cmdName));
	return it == cmdMap.end() ? 0 : it->second;
}

Cmd& Cmd::Engine::findCmd(const std::string& cmdName)
{
	Cmd* cmd = searchCmd(cmdName);
	if (cmd != 0)
		return *cmd;
	cmd = searchAlias(cmdName);
	if (cmd != 0)
		return *cmd;

	throw(std::string("Command not found: ") + cmdName);
}

Cmd& Cmd::Engine::nextCmd(std::istream& args)
{
	std::string cmdName;
	args >> cmdName;
	return findCmd(cmdName);
}

void Cmd::Engine::addAlias(const std::string & cmdAlias, const std::string & cmdName)
{
	cmdAliases[cmdAlias] = cmdName;
}

Cmd* Cmd::Engine::searchAlias(const std::string& cmdAlias)
{
	Cmd::AliasDict::iterator it = cmdAliases.find(to_lower(cmdAlias));
	if (it != cmdAliases.end())
	{
		return searchCmd(it->second);
	}
	return 0;
}

void Cmd::Engine::help(std::string& cmdName)
{
	Cmd* cmd = searchCmd(cmdName);
	if (cmd == 0)
	{
		cmd = searchAlias(cmdName);
	}

	if (cmd == 0)
	{
		help();
	}
	else
	{
		cmd->completeHelp();
	}
}

void Cmd::Engine::initHistPath()
{
	const char* home = getenv("HOME");
	std::ostringstream iHistPath;
	if (home != NULL)
		iHistPath << home << '/';
	iHistPath << '.' << idStr() << ".cmd";
	histPath = iHistPath.str();
}

void Cmd::Engine::help()
{
	std::ios::fmtflags flags(std::cout.flags()); // saving format flags
	std::cout << "---- List of all commands ----" << std::endl;
	for (Cmd::Map::iterator it = cmdMap.begin(); it != cmdMap.end(); ++it)
	{
		std::cout << std::setfill(' ') << std::setw(15) << std::left << it->first << " - ";
		it->second->help();
	}
	std::cout << "------------------------------" << std::endl;
	std::cout.flags(flags); // restoring format flags
}

void Cmd::help(void)
{
	std::cout << comment << std::endl;
}

void Cmd::completeHelp(void)
{
	std::cout << "Descr:\t" << comment << std::endl;
	std::cout << "Usage:\t" << usage() << std::endl;
}

} //cli
} //mdaq
