#ifndef CCFG_MODEL_COMMAND_COMMANDS_H_INCLUDED
#define CCFG_MODEL_COMMAND_COMMANDS_H_INCLUDED

enum CommandId {
	RenameConfigNodeCommand,
	SetConfigNodeValueCommand,
	SetConfigNodeUnitCommand,
	SetConfigNodeRangeCommand,
	SetConfigNodeIncludeUrlCommand,
	InsertConfigNodeCommand,
	OverrideConfigNodeCommand,
	DeleteConfigNodeCommand,
	IncludeConfigCommand,
	ImportConfigCommand,
	CommandCount
};

#endif // CCFG_MODEL_COMMAND_COMMANDS_H_INCLUDED
