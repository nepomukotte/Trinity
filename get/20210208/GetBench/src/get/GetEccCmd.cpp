/* =====================================================================================================================
 * GetEccCmd.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::GetEccCmd
 * class mdaq::GetEccCmd::Engine
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

#include "GetEccCmd.h"
#include "mdaq/DefaultPortNums.h"
#include "get/rc/SM.h"
#include "CCfg/Io/Url.h"
#include "CCfg/Io/Document.h"
#include "CCfg/Io/BufferStorage.h"
#include "CCfg/Io/XmlFormat.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
namespace ba = boost::algorithm;
#include <iostream>
#include <iterator>
using mdaq::utl::CmdException;

namespace get
{

/**
 * This command loads the register configuration into the ECC.
 */
struct SetISPM : public GetEccCmd
{
	SetISPM(GetEccCmd::Engine& engine) : GetEccCmd(engine, "ispm-set", "<module> <ISP1> <ISP2>", "Sets ISPM register") {}
	void exec(std::istream & cmdArgs)
	{
		std::string moduleName, nameISP1, nameISP2;

		cmdArgs >> moduleName;
		cmdArgs >> nameISP1;
		cmdArgs >> nameISP2;
		ecc()->asadSetISPM(moduleName, nameISP1, nameISP2);
	}
};

struct InjectCharge : public GetEccCmd
{
	InjectCharge(GetEccCmd::Engine& engine) : GetEccCmd(engine, "asad-pulser-trigger", "<voltage [mV]>", "Injects the given charge in mV") {}
	void exec(std::istream & cmdArgs)
	{
		long value;
		cmdArgs >> value;
		ecc()->triggerAsAdPulser(value);
	}
};

struct StartAsAdPeriodicPulser : public GetEccCmd
{
	StartAsAdPeriodicPulser(GetEccCmd::Engine& engine) : GetEccCmd(engine, "asad-pulser-start", "<voltage [mV]> <period [ms]>", "Starts periodic pulses of the AsAd voltage generator") {}
	void exec(std::istream & cmdArgs)
	{
		long voltage;
		size_t period;
		cmdArgs >> voltage >> period;
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		ecc()->configureAsAdPeriodicPulser(voltage, period);
		ecc()->startAsAdPeriodicPulser();
	}
};

struct StopAsAdPeriodicPulser : public GetEccCmd
{
	StopAsAdPeriodicPulser(GetEccCmd::Engine& engine) : GetEccCmd(engine, "asad-pulser-stop", "", "Stops periodic pulses of the AsAd voltage generator") {}
	void exec(std::istream &)
	{
		ecc()->stopAsAdPeriodicPulser();
	}
};

struct ListAllConfig : public GetEccCmd
{
	ListAllConfig(GetEccCmd::Engine& engine) : GetEccCmd(engine, "config-list", "[<config-subtype>]", "Returns list of all configuration IDs that are available to the ECC") {}
	void exec(std::istream & cmdArgs)
	{
		std::string subConfigType;
		cmdArgs >> subConfigType;
		if (cmdArgs.fail() or subConfigType.empty())
		{
			std::string configIds = ecc()->getConfigIds();

			std::cout << "Configuration identifiers:\n" << configIds << std::endl;
		}
		else
		{
			std::string subConfigIds = ecc()->getSubConfigIds(subConfigType);
			std::cout << "Configuration identifiers of type '" << subConfigType << "':\n" << subConfigIds << std::endl;
		}

	}
};

struct GetConfig : public GetEccCmd
{
	GetConfig(GetEccCmd::Engine& engine) : GetEccCmd(engine, "config-get", "<config-id>", "Retrieves configuration of given identifier.") {}
	void exec(std::istream & cmdArgs)
	{
		std::string configId;
		std::getline(cmdArgs, configId);
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		// Read configuration
		std::string contents = ecc()->getConfigSet(configId);
		LOG_INFO() << "Retrieved version with ID '" << configId << "':\n";
		std::cout << contents << std::endl;
	}
};

struct AddConfig : public GetEccCmd
{
	AddConfig(GetEccCmd::Engine& engine) : GetEccCmd(engine, "config-add", "<config-subtype> <file> [<sub-config-id>]", "Adds a configuration to the run control configuration manager.") {}
	void exec(std::istream & cmdArgs)
	{
		std::string subConfigType, configPath, subConfigId;
		cmdArgs >> subConfigType;
		cmdArgs >> configPath;
		cmdArgs >> subConfigId;

		// Read configuration
		CCfg::Io::Url configUrl(configPath);
		CCfg::Io::Document doc;
		doc.loadFromUrl(configUrl);
		CCfg::Io::BufferStorage buffer;
		CCfg::Io::XmlFormat xml;
		xml.format(doc.getConfig(), buffer, CCfg::Io::NoOpIncludeHandler::getInstance());
		std::string subConfig;
		buffer.readAll(subConfig);

		// Add config
		subConfigId = ecc()->addSubConfig(subConfigType, subConfig, subConfigId);
		LOG_INFO() << "Added version '" << subConfigId << "' of configuration '" << subConfigType << "'.";
	}
};

struct Describe : public GetEccCmd
{
	Describe(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-describe", "<config-id>", "Sets up hardware topology") {}
	void exec(std::istream & cmdArgs)
	{
		std::string configId;
		cmdArgs >> configId;
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		ecc()->describe(configId);
	}
};

struct Prepare : public GetEccCmd
{
	Prepare(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-prepare", "<config-id>", "Prepares inter-board communications, clock synchronizations, calibrations") {}
	void exec(std::istream & cmdArgs)
	{
		std::string configId;
		cmdArgs >> configId;
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		ecc()->prepare(configId);
	}
};

struct Reset : public GetEccCmd
{
	Reset(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-reset", "", "Resets ECC state machine, allowing to start configuration from scratch") {}
	void exec(std::istream & cmdArgs)
	{
		std::string configId;
		cmdArgs >> configId;
		ecc()->reset();
	}
};

struct Configure : public GetEccCmd
{

	Configure(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-configure", "<config-id> <data links>", "Loads physics parameters into the registers of the boards") {}
	void exec(std::istream & cmdArgs)
	{
		std::string configId, dataLinks;
		cmdArgs >> configId;
		std::getline(cmdArgs, dataLinks);
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		ecc()->configure(configId, dataLinks);
	}
};

struct Start : public GetEccCmd
{
	Start(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-start", "", "Starts the run, i.e. allows hardware to deliver data") {}
	void exec(std::istream &)
	{
		ecc()->start();
	}
};

struct Stop : public GetEccCmd
{

	Stop(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-stop", "", "Ends the run, preventing hardware from delivering data") {}
	void exec(std::istream &)
	{
		ecc()->stop();
	}
};

struct Pause : public GetEccCmd
{
	Pause(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-pause", "", "Pauses the run, i.e. forbids hardware to deliver data") {}
	void exec(std::istream &)
	{
		ecc()->pause();
	}
};

struct Resume : public GetEccCmd
{
	Resume(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-resume", "", "Resumes the run, i.e. allows again hardware to deliver data") {}
	void exec(std::istream &)
	{
		ecc()->resume();
	}
};

struct Breakup : public GetEccCmd
{

	Breakup(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-breakup", "", "Invalidates physics parameters setup in the boards") {}
	void exec(std::istream &)
	{
		ecc()->breakup();
	}
};

struct GetEccStatus : public GetEccCmd
{
	GetEccStatus(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-status", "", "Returns status (state, transition, error) of ECC state machine") {}
	void exec(std::istream &)
	{
		rc::Status status;
		ecc()->getStatus(status);
		std::cout << "State: " << rc::SM::State(status.s) << " Error: " << rc::SM::Error(status.e) << std::endl;
	}
};

struct GetDataLinks : public GetEccCmd
{
	GetDataLinks(GetEccCmd::Engine& engine) : GetEccCmd(engine, "sm-links", "", "Returns links between data senders and data routers") {}
	void exec(std::istream &)
	{
		std::cout << "Data links:\n" << ecc()->getDataLinks() << std::endl;
	}
};

struct PLLConfigure : public GetEccCmd
{
	PLLConfigure(GetEccCmd::Engine& engine) : GetEccCmd(engine, "pll-configure", "ckwFreq_MHz", "Configures PLL for chosen CKW frequency") {}
	void exec(std::istream & cmdArgs)
	{
		float value;
		cmdArgs >> value;
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		ecc()->loadPLLConfig(value);
	}
};

struct PLLWriteRegister : public GetEccCmd
{
	PLLWriteRegister(GetEccCmd::Engine& engine) : GetEccCmd(engine, "pll-write", "<value>", "Writes single PLL register") {}
	void exec(std::istream & cmdArgs)
	{
		uint32_t value;
		cmdArgs.unsetf(std::ios::dec);
		cmdArgs >> value;
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		LOG_DEBUG() << "writePLLRegister(" << std::hex << std::showbase << value << std::dec << std::noshowbase << ')';
		ecc()->writePLLRegister(value);
	}
};

struct PLLReset : public GetEccCmd
{
	PLLReset(GetEccCmd::Engine& engine) : GetEccCmd(engine, "pll-reset", "", "Sends RESET to PLL") {}
	void exec(std::istream &)
	{
		ecc()->resetPLL();
	}
};

struct MutantLinkSetup : public GetEccCmd
{
	MutantLinkSetup(GetEccCmd::Engine& engine) : GetEccCmd(engine, "cobo-mutant-link-set", "<0|1|2|3>", "Sets up mode of link between current CoBo node and Mutant (0=Normal, 1=Disconnected, 2=GMC, 3=WSCA)") {}
	void exec(std::istream & cmdArgs)
	{
		size_t modeValue;
		cmdArgs >> modeValue;
		if (cmdArgs.fail() or modeValue > 3)
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		cobo::MutantLinkMode mode = static_cast< cobo::MutantLinkMode >(modeValue);
		ecc()->setupCoBoMutantLink(mode);
	}
};

struct MutantCheckBram : public GetEccCmd
{
	MutantCheckBram(GetEccCmd::Engine& engine) : GetEccCmd(engine, "mutant-bram-check", "<value> <regCount> <firstOffset>", "Tests MUTANT_B BRAM by writing and reading <regCount> words") {}
	void exec(std::istream & cmdArgs)
	{
		::mdaq::hw::DeviceDescription dev;
		dev.name = "CARD_B_BRAM";
		dev.registerAccess = "MemBus16"; // "MemBus"
		dev.registerWidth = 2; // 4
		dev.baseAddress = 0x80F10000; // 0x200000E4

		uint64_t value = 0;
		uint64_t firstOffset = 0;
		uint64_t regCount;
		cmdArgs.unsetf(std::ios::dec);
		cmdArgs >> value >> regCount >> firstOffset;
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		LOG_INFO() << "Checking " << regCount << " word(s) with value 0x" << std::hex << value << " starting at offset 0x" << firstOffset << " from base address 0x" << dev.baseAddress;
		std::string result = ecc()->checkDevice(dev, value, firstOffset, regCount);
		LOG_INFO() << "Result of BRAM check: " << result;
	}
};

struct MutantAlignTDC : public GetEccCmd
{
	MutantAlignTDC(GetEccCmd::Engine& engine) : GetEccCmd(engine, "mutant-align", "<GMC|WSCA>", "Triggers alignment procedure of current Mutant node.") {}
	void exec(std::istream & cmdArgs)
	{
		std::string modeString = "GMC";
		cmdArgs >> modeString;
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		mutant::Clock clock = ("WSCA" != modeString)?mutant::GMC:mutant::WSCA;
		LOG_DEBUG() << "Triggering alignment procedure with clock " << clock << ", i.e. " << modeString;
		mutant::AlignmentResults results;
		bool verbose=true;
		try
		{
			// Trigger alignment procedure
			ecc()->mutantAlignmentServer()->alignDelayTDC(clock, verbose, results);
		}
		catch (const mutant::AlignmentAlarm & e)
		{
			throw CmdException(std::string("Error during alignment procedure: ") + e.reason);
		}
		catch (const Ice::LocalException & e)
		{
			throw CmdException(std::string("Could not trigger alignment procedure: ") + e.ice_name());
		}

		// Display results
		mutant::AlignmentResults::const_iterator resultIter;
		LOG_INFO() << results.size() << " result(s)";
		for (resultIter = results.begin(); resultIter != results.end(); ++resultIter)
		{
			LOG_INFO() << "Slot " << resultIter->first << ": "
					<< " isAmcPresent=" << (*resultIter).second.isAmcPresent
					<< " deltaDelay=" << (*resultIter).second.amcDeltaDelay
					<< " progDelay=" << (*resultIter).second.amcProgDelay;
		}
	}
};


struct MutantMultiShelveAlign : public GetEccCmd
{
	MutantMultiShelveAlign(GetEccCmd::Engine& engine) : GetEccCmd(engine, "mutant-align-multishelve", "<1|2|3> <GMC|WSCA>", "Multi shelve alignment procedure of current Mutant node.") {}
	void exec(std::istream & cmdArgs)
	{
		size_t step = 1;
		cmdArgs >> step;
		if (cmdArgs.fail())
		{
			throw CmdException("Invalid command parameters");
		}
		std::string modeString = "GMC";
		cmdArgs >> modeString;
		if (cmdArgs.fail())
		{
			throw CmdException("Invalid command parameters");
		}
		mutant::Clock clock = ("WSCA" != modeString)?mutant::GMC:mutant::WSCA;
		LOG_DEBUG() << "Triggering step " << step << "/3 of multi-shelve alignment procedure " << " with clock " << clock << ", i.e. " << modeString;
		bool verbose = true;
		mutant::MultiShelveAlignmentResults results;
		try {
			if (step == 1)
			{
				ecc()->mutantAlignmentServer()->multiShelveAlignDelayStep1(clock, verbose, results);
			}
			else if (step == 2)
			{
				ecc()->mutantAlignmentServer()->multiShelveAlignDelayStep2(clock, verbose, results);
			}
			else if (step == 3)
			{
				ecc()->mutantAlignmentServer()->multiShelveAlignDelayStep3(clock, verbose, results);
			}
		}
		catch (const mutant::AlignmentAlarm & e)
		{
			throw CmdException(std::string("Error during alignment procedure: ") + e.reason);
		}
		catch (const Ice::LocalException & e)
		{
			throw CmdException(std::string("Could not trigger alignment procedure: ") + e.ice_name());
		}

		// Display results
		if (not results.comment.empty())
		{
			LOG_INFO() << results.comment;
		}
	}
};

struct SwitchToCoBoSlowCtrlMode : public GetEccCmd
{
	SwitchToCoBoSlowCtrlMode(GetEccCmd::Engine& engine) : GetEccCmd(engine, "cobo-sc-slow", "", "Switches current CoBo node to slow control mode.") {}
	void exec(std::istream & cmdArgs)
	{
		bool enabled = true;
		cmdArgs >> std::boolalpha >> enabled;
		ecc()->coboSwitchToSlowControlMode(enabled);
	}
};

struct SwitchToCoBoFastCtrlMode : public GetEccCmd
{
	SwitchToCoBoFastCtrlMode(GetEccCmd::Engine& engine) : GetEccCmd(engine, "cobo-sc-fast", "", "Switches current CoBo node to hit register mode.") {}
	void exec(std::istream & cmdArgs)
	{
		bool enabled = true;
		cmdArgs >> std::boolalpha >> enabled;
		ecc()->coboSwitchToSlowControlMode(not enabled);
	}
};

struct SelectSingleAgetReadoutChannel : public GetEccCmd
{
	SelectSingleAgetReadoutChannel(GetEccCmd::Engine& engine) : GetEccCmd(engine, "aget-readout-channel-select", "<channel>", "Selects a single channel per chip to always read.") {}
	void exec(std::istream & cmdArgs)
	{
		size_t chanId = 68;
		cmdArgs >> chanId;
		ecc()->coboSelectSingleAgetReadoutChannel(chanId);
	}
};

struct SelectSingleAgetTestChannel : public GetEccCmd
{
	SelectSingleAgetTestChannel(GetEccCmd::Engine& engine) : GetEccCmd(engine, "aget-test-channel-select", "<channel>", "Selects a single channel per chip for test modes.") {}
	void exec(std::istream & cmdArgs)
	{
		size_t chanId = 68;
		cmdArgs >> chanId;
		ecc()->coboSelectSingleAgetTestChannel(chanId);
	}
};

struct SetAgetChannelLsbThreshold : public GetEccCmd
{
	SetAgetChannelLsbThreshold(GetEccCmd::Engine& engine) : GetEccCmd(engine, "aget-channel-lsb-treshold-set", "<channel> <threshold>", "Sets LSB threshold of a single channel per chip.") {}
	void exec(std::istream & cmdArgs)
	{
		size_t chanId = 68;
		cmdArgs >> chanId;
		uint16_t threshold = 15;
		cmdArgs >> threshold;
		ecc()->coboSetAgetChannelLsbThreshold(chanId, threshold);
	}
};

struct MonitorAsAdAlarms : public GetEccCmd
{
	MonitorAsAdAlarms(GetEccCmd::Engine& engine) : GetEccCmd(engine, "asad-alarm-monitor", "<true|false>", "Starts or stops the periodic check of the AsAd alarm bits.") {}
	void exec(std::istream & cmdArgs)
	{
		bool enabled = true;
		cmdArgs >> std::boolalpha >> enabled;
		ecc()->setAsAdAlarmMonitoringEnabled(enabled);
	}
};

struct SubscribeAlarms : public GetEccCmd
{
	SubscribeAlarms(GetEccCmd::Engine& engine) : GetEccCmd(engine, "node-alarm-subscribe", "<true|false>", "Subscribes ECC to alarms from current node.") {}
	void exec(std::istream & cmdArgs)
	{
		bool enabled = true;
		cmdArgs >> std::boolalpha >> enabled;
		ecc()->subscribeAlarms(enabled);
	}
};

struct BlinkLeds : public GetEccCmd
{
	BlinkLeds(GetEccCmd::Engine& engine) : GetEccCmd(engine, "cobo-led-blink-all", "", "Blinks LEDs.") {}
	void exec(std::istream &)
	{
		ecc()->ledManager()->blinkLEDs();
	}
};

struct SetLed : public GetEccCmd
{
	SetLed(GetEccCmd::Engine& engine) : GetEccCmd(engine, "cobo-led-set", "<ledId> <true|false>", "Switches LED on or off.") {}
	void exec(std::istream & cmdArgs)
	{
		size_t ledId;
		cmdArgs >> ledId;
		bool enabled = true;
		cmdArgs >> std::boolalpha >> enabled;
		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException(std::string("Invalid command parameters! Usage: ") + usage());
		}
		if (ledId == 16)
		{
			LOG_DEBUG() << "Switching " << (enabled?"on":"off") << " all LEDs";
			ecc()->ledManager()->setLEDs(enabled);
		}
		else
		{
			LOG_DEBUG() << "Switching " << (enabled?"on":"off") << " LED no. " << ledId;
			ecc()->ledManager()->setLED(ledId, enabled);
		}
	}
};

struct FlipLed : public GetEccCmd
{
	FlipLed(GetEccCmd::Engine& engine) : GetEccCmd(engine, "cobo-led-flip", "<ledId>", "Flips state of LED.") {}
	void exec(std::istream & cmdArgs)
	{
		size_t ledId;
		cmdArgs >> ledId;
		if (cmdArgs.fail() or ledId >= 16)
		{
			LOG_DEBUG() << "Flipping all LEDs";
			ecc()->ledManager()->flipLEDs();
		}
		else
		{
			LOG_DEBUG() << "Flipping LED no. " << ledId;
			ecc()->ledManager()->flipLED(ledId);
		}
	}
};

struct PulseLed : public GetEccCmd
{
	PulseLed(GetEccCmd::Engine& engine) : GetEccCmd(engine, "cobo-led-pulse", "<ledId> <period_ms>", "Flips LED periodically.") {}
	void exec(std::istream & cmdArgs)
	{
		size_t ledId, period_ms;
		cmdArgs >> ledId >> period_ms;
		if (cmdArgs.fail() or ledId > 16)
		{
			throw mdaq::utl::CmdException(std::string("Invalid command parameters! Usage: ") + usage());
		}
		else
		{
			if (ledId >= 16)
				LOG_DEBUG() << "Pulsing all LEDs";
			else
				LOG_DEBUG() << "Pulsing LED no. " << ledId;
			ecc()->ledManager()->pulseLED(ledId, period_ms);
		}
	}
};

GetEccCmd::Engine::Engine(const std::string & name, const std::string & servantName) : ::mdaq::EccCmd< GetEccPrx >::Engine(name, servantName)
{
	subscribe<SetISPM>();
	subscribe<InjectCharge>();
	subscribe<StartAsAdPeriodicPulser>();
	subscribe<StopAsAdPeriodicPulser>();
	subscribe<MonitorAsAdAlarms>();
	subscribe<SubscribeAlarms>();
	subscribe<ListAllConfig>();
	subscribe<GetConfig>();
	subscribe<AddConfig>();
	subscribe<Describe>();
	subscribe<Prepare>();
	subscribe<Reset>();
	subscribe<Configure>();
	subscribe<Start>();
	subscribe<Stop>();
	subscribe<Pause>();
	subscribe<Resume>();
	subscribe<Breakup>();
	subscribe<GetEccStatus>();
	subscribe<GetDataLinks>();
	subscribe<PLLConfigure>();
	subscribe<PLLWriteRegister>();
	subscribe<PLLReset>();
	subscribe<BlinkLeds>();
	subscribe<SetLed>();
	subscribe<FlipLed>();
	subscribe<PulseLed>();
	subscribe<MutantLinkSetup>();
	subscribe<MutantAlignTDC>();
	subscribe<MutantMultiShelveAlign>();
	subscribe<MutantCheckBram>();
	subscribe<SwitchToCoBoFastCtrlMode>();
	subscribe<SwitchToCoBoSlowCtrlMode>();
	subscribe<SelectSingleAgetReadoutChannel>();
	subscribe<SelectSingleAgetTestChannel>();
	subscribe<SetAgetChannelLsbThreshold>();

	addAlias("set-ispm", "ispm-set");
	addAlias("sm-topology", "sm-links");
	addAlias("asad-inject", "asad-pulser-trigger");
	addAlias("mutant-link-set", "cobo-mutant-link-set");
	addAlias("cobo-test-channel-select", "aget-test-channel-select");
	addAlias("cobo-readout-channel-select", "aget-readout-channel-select");
}

} // namespace get
