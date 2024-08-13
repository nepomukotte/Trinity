[Setup]
AppId={{361D7888-B4B7-11DE-C68D-0011439ED57D}
AppName=CompoundConfig
AppVerName=CompoundConfig 2.1.0
AppPublisher=Cea
AppPublisherURL=http://www.cea.fr
AppSupportURL=http://www.cea.fr
AppUpdatesURL=http://www.cea.fr
DefaultDirName={pf}\CompoundConfig
DefaultGroupName=CompoundConfig
OutputBaseFilename=setup
Compression=lzma
SolidCompression=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "french"; MessagesFile: "compiler:Languages\French.isl"

[Files]
Source: "_pkginst\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

