///
///
/// @brief A new file
///
/// @file xmlconst.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 26 mar 2022 
///

#ifndef XMLCONST_HPP_
#define XMLCONST_HPP_
#include <string>
#include <string_view>

namespace fcutils{
namespace monitor{

namespace xml {

	namespace cmd{
		typedef std::string key;
		namespace field {
			const key IDFLD="id";
			const key CMDFLD="cmd";
			const key HASDATAFLD="hasdata";
			const key CALCCRCFLD="calccrc";
			const key CMDSFLD="cmds";
		}
	}
	namespace seq{
		typedef std::string key;
		namespace field{
			const key IDFLD="id";
			const key SEQSFLD="seqs";
			const key SEQFLD="seq";
			const key CMDFLD="cmd";
			const key CMDSFLD="cmds";
		}
	}
	namespace client{
		typedef std::string key;
		namespace field{
			const key OPTSFLD="opts";
			const key OPTFLD="opt";
		}
		namespace option{
			const key IPOPT="clientServerIpAddress";
			const key HKCMDOPT="clientHkCommandFile";
			const key HKNCMDOPT="clientHkNCommandFile";
			const key NISTCMDOPT="clientNistCommandFile";
			const key PROCCMDOPT="clientProcCommandFile";
			const key SEQCMDOPT="clientSeqCommandFile";
			const key UCIRCCMDOPT="clientUcircCommandFile";
		}
	}

	namespace monitor{
		typedef std::string key;
		const key MONITORROOT="monitor";
		const key FAVORITESROOT="monitor.favorites";
		namespace field{
			const key INFOOPT="info";
			const key TITLEOPT="title";
			const key COUMNOPT="column";
			const key FAVORITEFILEOPT="file";
			const key FAVORITEOPT="favorite";
			const key FAVORITENAMEOPT="name";
		}
	}

	namespace slowcontrol{
		using key=std::string;
		using svkey=std::string_view;
		const key SLOWCONTROLROOT="slowcontrol";
		const key FAVORITESROOT="slowcontrol.favorites";
		const key VARIABLESROOT="slowcontrol.variables";
		namespace field{
			const key SUBSYSTEMOPT="subsystem";
			const key IDOPT="id";
			const key INFOOPT="info";
			const key PREFIXOPT="prefix";
			const key TEXTCOLOROPT="textcolor";
			const key BACKGROUNDCOLOROPT="backgroundcolor";
		}
		namespace value{
			const key SNONE="none";
			const key SRED="red";
			const key SGREEN="green";
			const key SBLUE="blue";
			const key SWHITE="white";
			const key SBLACK="black";
		}

		namespace variable{
			const key VARIABLEROOT="variables";
			const key VARIABLENODE="variable";
			namespace field{
				const key NAMEOPT="name";
				const key IDOPT="id";
				const key LABELOPT="label";
				const key OFFSETOPT="offset";
				const key SIZEOPT="size";
				const key INTYPEOPT="intype";
				const key OUTTYPEOPT="outtype";
				const key MAXOPT="max";
				const key MINOPT="min";
				const key INFOOPT="info";
				const key CONVERSIONFACTOROPT="conversion.factor";
				const key CONVERSIONOFFSETOPT="conversion.offset";
				const key CONVERSIONBYTEMASKOPT="conversion.bytemask";
				const key CONVERSIONBITSHIFTOPT="conversion.bitshift";
				const key ENDIANNESSOPT="endianness";
				const key FORMATOPT="format";
				const key ALARMCOLOROPT="alarmcolor";
			}
			namespace value{
				const key SBIGENDIAN="bigendian";
				const key SLITTLEENDIAN="littleendian";
				const key SDEC="dec";
				const key SHEX="hex";
				const key SINT="int";
				const key SFLOAT="float";
				const key SDOUBLE="double";
				const key SBOOL="bool";
				const key SCHAR="char";
				const key SSTRING="string";
				const key SNONE="none";
				const key SRED="red";
				const key SGREEN="green";
				const key SBLUE="blue";
				const key SWHITE="white";
				const key SBLACK="black";
				const key SDEFAULTKEY="def";
				const key DEFBYTEMASK="0xFFFFFFFFFFFFFFFF";
				constexpr const svkey SVBIGENDIAN="bigendian";
				constexpr const svkey SVLITTLEENDIAN="littleendian";
				constexpr const svkey SVDEC="dec";
				constexpr const svkey SVHEX="hex";
				constexpr const svkey SVINT="int";
				constexpr const svkey SVFLOAT="float";
				constexpr const svkey SVDOUBLE="double";
				constexpr const svkey SVBOOL="bool";
				constexpr const svkey SVCHAR="char";
				constexpr const svkey SVSTRING="string";
				constexpr const svkey SVNONE="none";
				constexpr const svkey SVRED="red";
				constexpr const svkey SVGREEN="green";
				constexpr const svkey SVBLUE="blue";
				constexpr const svkey SVWHITE="white";
				constexpr const svkey SVBLACK="black";
				constexpr const svkey SVDEFAULTKEY="def";
				const std::size_t DEFBITSHIFT=0;
			}
		}
		namespace favorite{
			const key FAVORITEROOT="favorites";
			const key FAVORITENODE="favorite";
			namespace field{
				const key NAMEOPT="name";
				const key INFOOPT="info";
				const key PRECISIONOPT="precision";
				const key WIDTHOPT="width";
				const key TEXTCOLOROPT="textcolor";
				const key BACKGROUNDCOLOROPT="backgroundcolor";
			}
			namespace value{
				const key DEFAULTNAME="SCDefault";
				const key SNONE="none";
				const key SRED="red";
				const key SGREEN="green";
				const key SBLUE="blue";
				const key SWHITE="white";
				const key SBLACK="black";
			}
		}
	} // namespace slowcontrol

	namespace gseserial{
		typedef std::string key;
		namespace field{
			const key OPTSFLD="opts";
			const key OPTFLD="opt";
		}
		namespace option{
			const key SERIALOPTION="gseserialSerialOptionFile";
			const key SERIALPORT="gseserialSerialPort";
			const key SERIALPORT2="gseserialSerialPort2";
			const key SHAREDMEMORYNAME="sharedMemoryName";
			const key PROCESSIDNAME="processID";
		}
	}
	namespace serialport{
		typedef std::string key;
		namespace option{
			const key OPTSFLD="ports";
			const key OPTFLD="port";
		}

		namespace field{
			const key NAMEOPT="name";
			const key BAUDOPT="baud";
			const key FLOWCTRLOPT="flowctrl";
			const key PARITYOPT="parity";
			const key STOPBITOPT="stopbit";
			const key CHARSIZEOPT="charsize";
		}
	}

	namespace envelope{
		typedef std::string key;
		namespace option{
			const key OPTSFLD="packets";
			const key OPTFLD="packet";
		}
		namespace packet{
			namespace field{
				const key IDOPT="id";
				const key HEADEROPT="header";
				const key TRAILEROPT="trailer";
				const key CRCOPT="crc";
				const key PADDINGOPT="padding";
			}
			namespace header{
				namespace field{
					const key SUBPACKETOPT="subpacket";
				}
				namespace subpacket{
					const key SIZEOPT="size";
					const key VALUEOPT="value";
				}
			}
			namespace trailer{
				namespace field{
					const key SUBPACKETOPT="subpacket";
				}
				namespace subpacket{
					const key SIZEOPT="size";
					const key VALUEOPT="value";
				}
			}
			namespace padding{
				namespace field{
					const key TYPEOPT="type";
					const key MAXBYTEOPT="max_byte";
				}
			}
		}
	}
}

} // namespace monitor
} // namespace fcutils



#endif // XMLCONST_HPP_ 
