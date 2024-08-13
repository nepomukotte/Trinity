module CCfg {
	sequence<string> RAttributeList;
	
	//CConfig Exceptions through Ice
	exception	ExceptionForIce
	{ string reason; };
	
	interface RCConfig {
		//GET FUNCTIONS
		string	getCurrentPath() throws ExceptionForIce;
		int		getInt(string item) throws ExceptionForIce;
		string	getStr(string item) throws ExceptionForIce;
		bool	getBool(string item) throws ExceptionForIce;
		float	getReal(string item) throws ExceptionForIce;
		string	getName(string item) throws ExceptionForIce;
		string	getIndex(string item) throws ExceptionForIce;
		string	getRange(string item) throws ExceptionForIce;
		string	getUnit(string item) throws ExceptionForIce;
		string	getPreComment(string item) throws ExceptionForIce;
		string	getPostComment(string item) throws ExceptionForIce;

		 //SET FUNCTIONS
		void	setCurrentPath(string path);
		void	setInt(string item, string value) throws ExceptionForIce;
		void	setStr(string item, string value) throws ExceptionForIce;
		void	setBool(string item, string value) throws ExceptionForIce;
		void	setReal(string item, string value) throws ExceptionForIce;
		void	setName(string item, string value) throws ExceptionForIce;
		void	setIndex(string item, string value) throws ExceptionForIce;
		void	setRange(string item, string value) throws ExceptionForIce;
		void	setUnit(string item, string value) throws ExceptionForIce;
		void	setPreComment(string item, string value) throws ExceptionForIce;
		void	setPostComment(string item, string value) throws ExceptionForIce;
 
		 //CONFIG FUNCTIONS
		void		append(string path) throws ExceptionForIce;
		RCConfig*	subConf(string path) throws ExceptionForIce;
		bool		testNode(string path) throws ExceptionForIce;
		void		createNode(string path) throws ExceptionForIce;
		void		removeNode(string path) throws ExceptionForIce;
		void		enumNode(string path, out RAttributeList node)throws ExceptionForIce;
		};
	 	
	 	//MANAGE DOCS
	 	interface	RDocument {
	 		RCConfig*	getConfig() throws ExceptionForIce;
	 		void		save()throws ExceptionForIce;
	 		void		saveToFile(string filePath)throws ExceptionForIce;
	 		void		saveToUrl(string url)throws ExceptionForIce;
	 	};
	 	
	 	//CREATE DOCS
	 	interface	ConfigServer {
	 		RDocument*	loadDoc(string fileName);
	 		RDocument*	newDoc();
	 	};
};
