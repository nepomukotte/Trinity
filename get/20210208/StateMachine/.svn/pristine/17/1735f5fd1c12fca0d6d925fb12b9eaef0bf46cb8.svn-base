module dhsm
{
	exception ServantError
	{
		string reason;
	};

	exception ServantCriticalError extends ServantError
	{
	};

	interface IEngineServant
	{
		["ami"] int start() throws ServantError;
		["ami"] int reinitiate(string stateMachine) throws ServantError;
		["ami"] int halt(string stateMachine) throws ServantError;

		["ami"] int stepEvent(string event) throws ServantError;
		["ami"] int stepStringMsg(string event, string param) throws ServantError;
		["ami"] int stepIntMsg(string event, int param) throws ServantError;
	};
};
