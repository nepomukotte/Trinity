/*
 * $Id: RegisterProcessor.h 1780 2016-05-03 15:31:16Z anvar $
 * @file RegisterProcessor.h
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

#ifndef mdaq_hw_RegisterProcessor_h_INCLUDED
#define mdaq_hw_RegisterProcessor_h_INCLUDED

//#include "mdaq/hw/RegCmd.h"

#include "mdaq/hw/Control.h"
#include "mdaq/utl/VersaPointer.hpp"
#include "utl/Logging.h"

#include ICEUTIL_HEADER(Time.h)
#include ICEUTIL_HEADER(Thread.h)

#include <boost/ptr_container/ptr_list.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/noncopyable.hpp>
#include <string>
#include <vector>
#include <map>

namespace mdaq
{
namespace hw
{

void initRegCmd(RegCmd& regCmd, RegInstruction instruction, const std::string& regName = "", const std::string& fieldName = "", long value = -1);
RegCmd makeRegCmd(RegInstruction instruction, const std::string& regName = "", const std::string& fieldName = "", long value = -1);
RegCmd* newRegCmd(RegInstruction instruction, const std::string& regName = "", const std::string& fieldName = "", long value = -1);

template< RegInstruction I> struct Cmd;

/**
 * This is a class representing a "register access sequence processor". It is passed a "core object" that implements the
 * nested Core interface. The Core class operations (beginning with prefix "impl_") are the elementary instructions
 * executed by any sequence of register read/write operations.
 */
class RegisterProcessor : boost::noncopyable
{
public:
	// Nested types forward declarations
	class Cursor;
	class Core;
	class DefaultCore;
	typedef std::auto_ptr<Core> CorePtr;
	class CmdExec;
	class BatchCore;
	typedef boost::ptr_map<std::string, mdaq::hw::RegisterProcessor> Map;

public:
	RegisterProcessor(CorePtr ap_core_= CorePtr());
	virtual ~RegisterProcessor();

	void setCore(CorePtr ap_core_);

	Core& core() { return *ap_core; }
	void init();
	void flush();

	bool isFlushed() const { return flushed; }

	template < RegInstruction I >
	RegisterProcessor& exec(
		const std::string& regName,
			const std::string& fieldName,
			Ice::Long inValue,
			Ice::Long& outValue
			);

	template < RegInstruction I >
	RegisterProcessor& exec(
		const std::string& regName,
			const std::string& fieldName,
			Ice::Long inValue
			);

	RegisterProcessor& exec(
			RegInstruction instr,
			const std::string& regName,
			const std::string& fieldName,
			Ice::Long inValue,
			Ice::Long& outValue
			);

	RegisterProcessor& exec(
			RegInstruction instr,
			const std::string& regName,
			const std::string& fieldName,
			Ice::Long inValue
			);

	RegisterProcessor& exec(const RegCmd& cmd, RegCmd& trace);
	RegisterProcessor& exec(const RegCmdSeq& cmdSeq, RegCmdSeq& trace);

	RegisterProcessor& operator << (void (RegisterProcessor::*func)())
	{
		if (func == &RegisterProcessor::init)
			init();
		else if (func == &RegisterProcessor::flush)
			flush();
		return *this;
	}
	RegisterProcessor& operator << (const Cmd< Nop >& cmd_);
	RegisterProcessor& operator << (const Cmd<Read >& cmd_);
	RegisterProcessor& operator << (const Cmd<Write>& cmd_);
	RegisterProcessor& operator << (const Cmd<Sleep>& cmd_);

	Cursor write(const std::string& registerName, const std::string& bitFieldName, Ice::Long value);
	Cursor write(const std::string& registerName, Ice::Long value);
	Cursor read(const std::string& registerName, const std::string& bitFieldName, Ice::Long& value);
	Cursor read(const std::string& registerName, Ice::Long& value);
	Cursor sleep_ms(Ice::Long value);
	Cursor nop();

private:
	CorePtr ap_core; ///< When Core object's life cycle should belong to RegisterProcessor.
	bool flushed; ///< This is to ensure that flush operation is executed only once per sequence.
};

/**
 * This class shares the same processing interface as RegisterProcessor and is used to chain successive calls to the
 * RegisterProcessor instance to which it is associated. The Cursor object is instantiated on the stack (and its
 * constructor called) by the call of the first instruction of the sequence on a RegisterProcessor object. The end
 * of the processing chain causes the Cursor destructor to be called, which in turn calls the
 * RegisterProcessor::flush virtual method. The implementation of the impl_init() and impl_flush methods of
 * RegisterProcessor::Core can be used to respectively initialize and finalize the whole sequence processing.
 *
 * Every instruction sequence begins with a call to RegisterProcessor::init (in Cursor constructor) and ends with a
 * call to RegisterProcessor::flush (in Cursor destructor).
 */
class RegisterProcessor::Cursor
{
public:
	Cursor(RegisterProcessor& registerProcessor_) : registerProcessor(registerProcessor_) { registerProcessor.init(); }
	~Cursor() { registerProcessor.flush(); }

	Cursor& write(const std::string& registerName, const std::string& bitFieldName, Ice::Long value);
	Cursor& write(const std::string& registerName, Ice::Long value);
	Cursor& read(const std::string& registerName, const std::string& bitFieldName, Ice::Long& value);
	Cursor& read(const std::string& registerName, Ice::Long& value);
	Cursor& sleep_ms(Ice::Long value);
	Cursor& nop();

private:
	RegisterProcessor& registerProcessor;
};

/**
 * This class is the abstract public interface that specialized classes must implement as their execution model.
 * The virtual functions that must be implemented execute the basic operations of a register processor.
 */
class RegisterProcessor::Core
{
public:
	Core() {}
	virtual ~Core() {}
	virtual void impl_init() {};
	virtual void impl_w(const std::string& registerName, const std::string& bitFieldName, Ice::Long value) = 0;
	virtual void impl_r(const std::string& registerName, const std::string& bitFieldName, Ice::Long& value) = 0;
	virtual void impl_sleep_ms(Ice::Long value) { IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(value)); }
	virtual void impl_nop() {};
	virtual void impl_flush() {};
};

/**
 * This class is a basic adapter of the Core the abstract public interface. It defines a basic register structure and
 * actually implements read/write methods that create specified registers over a map.
 */
class RegisterProcessor::DefaultCore : public Core
{
public:
	typedef std::map< std::string, std::map<std::string, long> > RegisterMap;

public:
	virtual void impl_w(const std::string& registerName, const std::string& bitFieldName, Ice::Long value);
	virtual void impl_r(const std::string& registerName, const std::string& bitFieldName, Ice::Long& value);

private:
	RegisterMap regMap;
};

class RegisterProcessor::CmdExec
{
public:
	CmdExec(RegInstruction instr = Nop, const std::string& regName_ = "", const std::string& fieldName_ = "",
			Ice::Long value_ = -1, Ice::Long* p_outValue_ = 0)
		: regCmd_(newRegCmd(instr, regName_, fieldName_, value_)), p_outValue(p_outValue_)
	{
		if (p_outValue == 0)
			p_outValue = &temp;
	}

	CmdExec(RegCmd& regCmd_, Ice::Long* p_outValue_ = 0) : regCmd_(regCmd_), p_outValue(p_outValue_)
	{
		if (p_outValue == 0)
			p_outValue = &temp;
	}

	virtual ~CmdExec() {}

	virtual void exec(RegisterProcessor::Core& core) = 0;

	const RegCmd   & regCmd  () const { return *regCmd_   ; }
	      RegCmd   & regCmd  ()       { return *regCmd_   ; }

	      Ice::Long& outValue() const { return *p_outValue; }

private:
	mutable mdaq::utl::VersaPointer<RegCmd> regCmd_;
	mutable Ice::Long* p_outValue;
	Ice::Long temp;
};

template <>
struct Cmd<Nop> : public RegisterProcessor::CmdExec
{
	Cmd() : CmdExec() {}
	Cmd(RegCmd& regCmd_) : CmdExec(regCmd_) {}
	void exec(RegisterProcessor::Core& core) { core.impl_nop(); }
};

template <>
struct Cmd<Read> : public RegisterProcessor::CmdExec
{
	Cmd(const std::string& regName, const std::string& fieldName, Ice::Long* p_outValue_ = 0)
			: RegisterProcessor::CmdExec(Read, regName, fieldName, -1, p_outValue_) {}

	Cmd(RegCmd& regCmd_, Ice::Long* p_outValue_) : RegisterProcessor::CmdExec(regCmd_, p_outValue_) {}

	void exec(RegisterProcessor::Core& core) { core.impl_r(regCmd().regName, regCmd().fieldName, outValue()); }
};

template <>
struct Cmd<Write> : public RegisterProcessor::CmdExec
{
	Cmd(const std::string& regName, const std::string& fieldName, Ice::Long inValue)
			: RegisterProcessor::CmdExec(Write, regName, fieldName, inValue) {}

	Cmd(RegCmd& regCmd_) : RegisterProcessor::CmdExec(regCmd_) {}

	void exec(RegisterProcessor::Core& core) { core.impl_w(regCmd().regName, regCmd().fieldName, regCmd().value); }
};

template <>
struct Cmd<Sleep> : public RegisterProcessor::CmdExec
{
	Cmd(Ice::Long value_ms) : RegisterProcessor::CmdExec(Sleep, "", "", value_ms) {}
	Cmd(RegCmd& regCmd_) : RegisterProcessor::CmdExec(regCmd_) {}
	void exec(RegisterProcessor::Core& core) { core.impl_sleep_ms(regCmd().value); }
};

/**
 * This class is a specialization of Core that is used to substitute the Core object of a register processor by
 * a batched version of it, i.e. turning the execution model into a batched execution model: the whole sequence is
 * accumulated inside a command sequence that is eventually executed in the flush method.
 */
class RegisterProcessor::BatchCore : public Core
{
public:
	typedef boost::ptr_list<CmdExec> CmdSeq;

public:
	BatchCore(CorePtr ap_core = CorePtr());
	void impl_init();
	void impl_w(const std::string& regName, const std::string& fieldName, Ice::Long value);
	void impl_r(const std::string& regName, const std::string& fieldName, Ice::Long& value);
	void impl_sleep_ms(Ice::Long value);
	void impl_nop();
	void impl_flush();

	Core& actualCore() { return *actualCore_; }
	void getRegCmdSeq(RegCmdSeq& regCmdSeq) const;
	RegCmdSeq regCmdSeq() const;

protected:
	CmdSeq cmdBatch;
	CorePtr actualCore_;
};

} /* namespace hw */
} /* namespace mdaq */

#endif /* mdaq_hw_RegisterProcessor_h_INCLUDED */
