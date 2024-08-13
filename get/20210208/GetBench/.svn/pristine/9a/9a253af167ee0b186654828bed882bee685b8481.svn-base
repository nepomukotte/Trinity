//gsoap ecc  service name:	ecc 
//gsoap ecc  service port:	http://localhost:8061 
//gsoap ecc  service namespace:	urn:ecc
//gsoap ecc  service encoding: literal
//gsoap ecc  service transport:	http://schemas.xmlsoap.org/soap/http 
//gsoap ecc  service style:    document

class ecc__Response
{
  int         ErrorCode;
  std::string ErrorMessage;
};

class ecc__ResponseState
{
  int         ErrorCode;
  std::string ErrorMessage;
  int         State;
  int         Transition;
};

class ecc__ResponseText
{
  int         ErrorCode;
  std::string ErrorMessage;
  std::string Text;
};

int ecc__ServerExit (ecc__Response & response);

int ecc__GetConfigIDs (ecc__ResponseText & response);

int ecc__Describe (std::string configID, std::string table, ecc__ResponseText & response);

int ecc__Prepare (std::string configID, std::string table, ecc__ResponseText & response);

int ecc__Configure (std::string configID, std::string table, ecc__ResponseText & response);

int ecc__Start (std::string configID, std::string table, ecc__ResponseText & response);

int ecc__Stop (std::string configID, std::string table, ecc__ResponseText & response);

int ecc__Pause (std::string configID, std::string table, ecc__ResponseText & response);

int ecc__Resume (std::string configID, std::string table, ecc__ResponseText & response);

int ecc__Breakup (std::string configID, std::string table, ecc__ResponseText & response);

int ecc__Undo (std::string configID, std::string table, ecc__ResponseText & response);

int ecc__GetState (ecc__ResponseState & response);

int ecc__Special (std::string command, std::string table, ecc__ResponseText & response);
