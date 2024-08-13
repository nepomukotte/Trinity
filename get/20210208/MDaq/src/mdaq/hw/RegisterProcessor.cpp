/*
 * $Id: RegisterProcessor.cpp 1780 2016-05-03 15:31:16Z anvar $
 * @file RegisterProcessor.cpp
 * @created May 3, 2016
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of the Mordicus Real-Time Software
 * Development Framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free
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

#include "mdaq/hw/RegisterProcessor.h"
#include "mdaq/hw/Control.h"
#include <utl/Logging.h>
#include <time.h>

namespace mdaq
{
namespace hw
{

/**
 * This is a helper function that makes up for the absence of a constructor in Ice generated structure RegCmd.
 * @param regCmd The RegCmd instance to be initialized.
 * @param instruction The instruction that this command must execute
 * @param regName The name of the register to which this command is applied (when needed).
 * @param fieldName The name of the register field to which this command is applied (when needed).
 * @param value The value parameter associated to this command's instruction.
 */
void initRegCmd(RegCmd& regCmd, RegInstruction instruction, const std::string& regName, const std::string& fieldName, long value)
{
	regCmd.instruction = instruction;
	regCmd.regName = regName;
	regCmd.fieldName = fieldName;
	regCmd.value = value;
}

/**
 * This is a helper function that makes up for the absence of a constructor in Ice generated structure RegCmd.
 * @param instruction The instruction that this command must execute
 * @param regName The name of the register to which this command is applied (when needed).
 * @param fieldName The name of the register field to which this command is applied (when needed).
 * @param value The value parameter associated to this command's instruction.
 * @return By value, a properly initialized RegCmd copy.
 */
RegCmd makeRegCmd(RegInstruction instruction, const std::string& regName, const std::string& fieldName, long value)
{
	RegCmd regCmd;
	initRegCmd(regCmd, instruction, regName, fieldName, value);
	return regCmd;
}

/**
 * This is a helper function that makes up for the absence of a constructor in Ice generated structure RegCmd.
 * @param instruction The instruction that this command must execute
 * @param regName The name of the register to which this command is applied (when needed).
 * @param fieldName The name of the register field to which this command is applied (when needed).
 * @param value The value parameter associated to this command's instruction.
 * @return A newly allocated, properly initialized RegCmd (ownership is left to caller).
 */
RegCmd* newRegCmd(RegInstruction instruction, const std::string& regName, const std::string& fieldName, long value)
{
	RegCmd* p_regCmd = new RegCmd;
	initRegCmd(*p_regCmd, instruction, regName, fieldName, value);
	return p_regCmd;
}

/**
 * Constructs an RegisterProcessor in charge of executing register access operations.
 * @param core Auto pointer to the core object implementing the basic operations of a register processor; this object
 *        implements all the operations of the RegisterProcessor::Core abstract interface. The auto pointer is released
 *        and the register processor becomes the owner of the allocated Core object.
 */
RegisterProcessor::RegisterProcessor(CorePtr ap_core_) : flushed(true)
{
	setCore(ap_core_);
}

/**
 * If RegisterProcessor hasn't been flushed, flush() is called.
 */
RegisterProcessor::~RegisterProcessor()
{
	if (not flushed)
	{
		flush();
	}
}

void RegisterProcessor::setCore(CorePtr ap_core_)
{
	ap_core.reset(ap_core_.release());
	if (ap_core.get() == 0)
	{
		ap_core.reset(new DefaultCore);
	}
}

void RegisterProcessor::init()
{
	core().impl_init();
	flushed = false;
}

void RegisterProcessor::flush()
{
	flushed = true; // Must be set before actually calling flush in case execution of flush causes exception.
	core().impl_flush();
}

template <>
RegisterProcessor&
RegisterProcessor::exec<Nop> (const std::string&, const std::string&, Ice::Long)
{
	core().impl_nop();
	return *this;
}

template <>
RegisterProcessor&
RegisterProcessor::exec<Write> (const std::string& regName, const std::string& fieldName, Ice::Long inValue)
{
	core().impl_w(regName, fieldName, inValue);
	return *this;
}

template <>
RegisterProcessor&
RegisterProcessor::exec<Read> (const std::string& regName, const std::string& fieldName, Ice::Long, Ice::Long& outValue)
{
	core().impl_r(regName, fieldName, const_cast<Ice::Long&>(outValue));
	return *this;
}

template <>
RegisterProcessor&
RegisterProcessor::exec<Sleep> (const std::string&, const std::string&, Ice::Long value_ms)
{
	core().impl_sleep_ms(value_ms);
	return *this;
}

RegisterProcessor&
RegisterProcessor::exec (
		RegInstruction instruction,
		const std::string& regName,
		const std::string& fieldName,
		Ice::Long inValue,
		Ice::Long& outValue
		)
{
	switch(instruction)
	{
	case Write:
		exec<Write>(regName, fieldName, inValue);
		break;

	case Read:
		exec<Read>(regName, fieldName, inValue, outValue);
		break;

	case Sleep:
		exec<Sleep>(regName, fieldName, inValue);
		break;

	case Nop:
		exec<Nop>(regName, fieldName, inValue);
		break;
	}
	return *this;
}

RegisterProcessor& RegisterProcessor::exec(const RegCmd& cmd, RegCmd& trace)
{
	if (&trace != &cmd) trace = cmd;
	return exec(cmd.instruction, cmd.regName, cmd.fieldName, cmd.value, trace.value);
}

RegisterProcessor& RegisterProcessor::exec(const RegCmdSeq& cmdSeq, RegCmdSeq& trace)
{
	trace.clear();
	for(RegCmdSeq::const_iterator it = cmdSeq.begin(); it != cmdSeq.end(); ++it)
	{
		trace.push_back(*it);
		exec(*it, trace.back());
	}
	return *this;
}

RegisterProcessor& RegisterProcessor::operator << (const Cmd< Nop >& cmd_)
{
	exec< Nop >(cmd_.regCmd().regName, cmd_.regCmd().fieldName, cmd_.regCmd().value);
	return *this;
}

RegisterProcessor& RegisterProcessor::operator << (const Cmd<Read >& cmd_)
{
	exec<Read >(cmd_.regCmd().regName, cmd_.regCmd().fieldName, cmd_.regCmd().value, cmd_.outValue());
	return *this;
};

RegisterProcessor& RegisterProcessor::operator << (const Cmd<Write>& cmd_)
{
	exec<Write>(cmd_.regCmd().regName, cmd_.regCmd().fieldName, cmd_.regCmd().value);
	return *this;
}

RegisterProcessor& RegisterProcessor::operator << (const Cmd<Sleep>& cmd_)
{
	exec<Sleep>(cmd_.regCmd().regName, cmd_.regCmd().fieldName, cmd_.regCmd().value);
	return *this;
}

RegisterProcessor::Cursor
RegisterProcessor::write(const std::string& regName, const std::string& fieldName, Ice::Long value)
{
	Cursor cursor(*this);
	exec<Write>(regName, fieldName, value);
	return cursor;
}

RegisterProcessor::Cursor
RegisterProcessor::write(const std::string& regName, Ice::Long value)
{
	Cursor cursor(*this);
	exec<Write>(regName, "", value);
	return cursor;
}

RegisterProcessor::Cursor
RegisterProcessor::read(const std::string& regName, const std::string& fieldName, Ice::Long& value)
{
	Cursor cursor(*this);
	exec<Read>(regName, fieldName, value, value);
	return cursor;
}

RegisterProcessor::Cursor
RegisterProcessor::read(const std::string& regName, Ice::Long& value)
{
	Cursor cursor(*this);
	exec<Read>(regName, "", value, value);
	return cursor;
}

RegisterProcessor::Cursor
RegisterProcessor::sleep_ms(Ice::Long value)
{
	Cursor cursor(*this);
	exec<Sleep>("", "", value);
	return cursor;
}

RegisterProcessor::Cursor
RegisterProcessor::nop()
{
	Ice::Long fakeValue = 0;
	Cursor cursor(*this);
	exec<Nop>("", "", fakeValue);
	return cursor;
}

RegisterProcessor::Cursor&
RegisterProcessor::Cursor::write(const std::string& regName, const std::string& fieldName, Ice::Long value)
{
	registerProcessor.exec<Write>(regName, fieldName, value);
	return *this;
}

RegisterProcessor::Cursor&
RegisterProcessor::Cursor::write(const std::string& regName, Ice::Long value)
{
	registerProcessor.exec<Write>(regName, "", value);
	return *this;
}

RegisterProcessor::Cursor&
RegisterProcessor::Cursor::read(const std::string& regName, const std::string& fieldName, Ice::Long& value)
{
	registerProcessor.exec<Read>(regName, fieldName, value, value);
	return *this;
}

RegisterProcessor::Cursor&
RegisterProcessor::Cursor::read(const std::string& regName, Ice::Long& value)
{
	registerProcessor.exec<Read>(regName, "", value, value);
	return *this;
}

RegisterProcessor::Cursor&
RegisterProcessor::Cursor::sleep_ms(Ice::Long value)
{
	registerProcessor.exec<Sleep>("", "", value);
	return *this;
}

RegisterProcessor::Cursor&
RegisterProcessor::Cursor::nop()
{
	Ice::Long fakeValue = 0;
	registerProcessor.exec<Nop>("", "", fakeValue);
	return *this;
}

void RegisterProcessor::DefaultCore::impl_w(const std::string& regName, const std::string& fieldName, Ice::Long value)
{
		regMap[regName][fieldName] = value;
}

void RegisterProcessor::DefaultCore::impl_r(const std::string& regName, const std::string& fieldName, Ice::Long& value)
{
		value = regMap[regName][fieldName];
}

RegisterProcessor::BatchCore::BatchCore(CorePtr ap_core)
{
	actualCore_.reset(ap_core.release());
	if (actualCore_.get() == 0)
		actualCore_.reset(new DefaultCore);
}

void RegisterProcessor::BatchCore::impl_init()
{
	cmdBatch.clear();
}

void RegisterProcessor::BatchCore::impl_w(const std::string& regName, const std::string& fieldName, Ice::Long value)
{
	cmdBatch.push_back(std::auto_ptr< CmdExec >(new Cmd<Write>(regName, fieldName, value)));
}

void RegisterProcessor::BatchCore::impl_r(const std::string& regName, const std::string& fieldName, Ice::Long& value)
{
	cmdBatch.push_back(std::auto_ptr< CmdExec >(new Cmd<Read>(regName, fieldName, &value)));
}

void RegisterProcessor::BatchCore::impl_sleep_ms(Ice::Long value_ms)
{
	cmdBatch.push_back(std::auto_ptr< CmdExec >(new Cmd<Sleep>(value_ms)));
}

void RegisterProcessor::BatchCore::impl_nop()
{
	cmdBatch.push_back(std::auto_ptr< CmdExec >(new Cmd<Nop>()));
}

void RegisterProcessor::BatchCore::impl_flush()
{
	if (actualCore_.get())
	{
		for (CmdSeq::iterator it = cmdBatch.begin(); it != cmdBatch.end(); ++it)
		{
			it->exec(actualCore());
		}
	}
}

void RegisterProcessor::BatchCore::getRegCmdSeq(RegCmdSeq& regCmdSeq) const
{
	regCmdSeq.resize(cmdBatch.size());
	CmdSeq::const_iterator itCmd = cmdBatch.begin();
	for (size_t i=0; i < regCmdSeq.size(); ++i)
	{
		regCmdSeq[i] = itCmd->regCmd();
		++itCmd;
	}
}

RegCmdSeq RegisterProcessor::BatchCore::regCmdSeq() const
{
	RegCmdSeq regCmdSeq;
	getRegCmdSeq(regCmdSeq);
	return regCmdSeq;
}

} /* namespace hw */
} /* namespace mdaq */
