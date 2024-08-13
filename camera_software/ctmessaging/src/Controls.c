#include "CTM_ProcessIDs.h"
#include "CTM_RemoteControlCommands.h"
#include "CTM_Controls.h"
#include "CTM_ConstantValues.h"
#include "CTM_TriggerBoardCommands.h"
#include "CTM_SIABCommands.h"
#include "CTM_CoBoCommands.h"
#include "CTM_PDUCommands.h"
#include "CTM_LVPSCommands.h"
#include "CTM_SequenceCommands.h"
#include "CTM_GyroCommands.h"
#include <mqueue.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <curl/curl.h>

//Navigation Routines
void Backspace(){
  char del[] = "\b \b";
  write(STDOUT_FILENO, del, sizeof(del) - 1);
}
void Clear(){
  char blank = malloc(STREAM_LENGTH);
  write(STDOUT_FILENO,blank,sizeof(blank)-1);
  //free(blank);
}
void Left(int x){
  printf("\033[%dD",x);
  fflush(stdout);
}
void Right(int x){
  printf("\033[%dC",x);
  fflush(stdout);
}


//Formatting
void PrintRed(){
  printf("%s",KRED);
}
void PrintBlue(){
  printf("%s",KBLU);
}
void PrintCyan(){
  printf("%s",KCYN);
}
void PrintGreen(){
  printf("%s",KGRN);
}
void PrintYellow(){
  printf("%s",KYLW);
}
void PrintMagenta(){
  printf("%s",KMAG);
}
void PrintWhite(){
  printf("%s",KWHT);
}
void PrintGray(){
  printf("%s",KGRY);
}
void PrintDefault(){
  printf("%s",KDEF);
}



void GetDictionaryEntry(char *dictEntry,const char* dict){
  sprintf(dictEntry, "%s", dict);
}

//String Routines
void RemoveTrailingCharacters(char* str){
  str[strcspn(str,"\r\n\b")]='\0';
}

void ClearString(char *str){
  memset(&str[0], 0, strlen(str)+1);
}

void CleanString(char *in){
  int i = 0;
  while(isalnum(in[i]) || in[i]==' ' || in[i] == '-'){
    i = i+1;
  }
  *(in+i) = '\0';

}

int FindSmaller(int num1, int num2){
  int r = 0;
  if(num1<num2){
    r=num1;
  }else{
    r=num2;
  }

  return r;
}

int SearchProcessID(char *pID){
  int i = 0;
  int isFound = 0;
  char strName[20];
  int comparisonSize=0;
  int pIDVal=-1;
  

  RemoveTrailingCharacters(pID);
  while(!isFound && i<NTOTALPROCESS){
    
    ClearString(strName);
    GetDictionaryEntry(strName,processID.names[i]); 
    sprintf(strName,"%s",processID.names[i]);   
    RemoveTrailingCharacters(strName);
    comparisonSize = FindSmaller(strlen(strName),strlen(pID));

    if(strncmp(pID,strName,comparisonSize) == 0){
      isFound = 1;
      pIDVal = processID.values[i];
      //printf("found target %d", i);
    }
    else{
      i = i+1;
    }
  };


  return pIDVal;
}

int FindRemoteControlCommand(char *cmd){
  int cmd_val;
  int isFound = 0;
  int i = 0;
  int comparisonSize = 0;
  char strName[STREAM_LENGTH];
  while(!isFound && i<NREMOTECTRLCMDS){
    ClearString(strName);
    GetDictionaryEntry(strName,remoteControlCommands.names[i]);    
    RemoveTrailingCharacters(strName);
    comparisonSize = FindSmaller(strlen(cmd),strlen(strName));
    if(strncmp(cmd,strName,comparisonSize)==0){
     isFound = 1;
     //printf("Here");
    }else{
     i = i +1;
    }
  };
  if(isFound){
    cmd_val = remoteControlCommands.values[i];
  }else{
    cmd_val = -1;
  }
  
  return cmd_val;
}

int FindTriggerBoardCommand(char *cmd){
  int cmd_val;
  int isFound = 0;
  int i = 0;
  int comparisonSize = 0;
  char strName[STREAM_LENGTH];
  while(!isFound && i<NTBCOMMANDS){
    ClearString(strName);
    GetDictionaryEntry(strName,triggerBoardCommands.names[i]);    
    RemoveTrailingCharacters(strName);
    //printf("%s %s\n",strName,cmd);
    //fflush(stdout);
    //comparisonSize = FindSmaller(strlen(cmd),strlen(strName));
    if(strcmp(cmd,strName)==0){
     isFound = 1;
     //printf("Here");
    }else{
     i = i +1;
    }
  };
  if(isFound){
    cmd_val = triggerBoardCommands.values[i];
  }else{
    cmd_val = -1;
  }
  
  return cmd_val;
}

int FindSIABCommand(char *cmd){
  int cmd_val;
  int isFound = 0;
  int i = 0;
  int comparisonSize = 0;
  char strName[STREAM_LENGTH];
  while(!isFound && i<NSIABCOMMANDS){
    ClearString(strName);
    GetDictionaryEntry(strName,siabCommands.names[i]);    
    RemoveTrailingCharacters(strName);
    comparisonSize = FindSmaller(strlen(cmd),strlen(strName));
    if(strncmp(cmd,strName,comparisonSize)==0){
     isFound = 1;
     //printf("Here");
    }else{
     i = i +1;
    }
  };
  if(isFound){
    cmd_val = siabCommands.values[i];
  }else{
    cmd_val = -1;
  }
  
  return cmd_val;
}

int FindCoBoCommand(char *cmd){
  int cmd_val;
  int isFound = 0;
  int i = 0;
  int comparisonSize = 0;
  char strName[STREAM_LENGTH];
  while(!isFound && i<NCOBOCOMMANDS){
    ClearString(strName);
    GetDictionaryEntry(strName,coboCommands.names[i]);    
    RemoveTrailingCharacters(strName);
    comparisonSize = FindSmaller(strlen(cmd),strlen(strName));
    if(strcmp(cmd,strName)==0){
     isFound = 1;
     //printf("%s %sHere",cmd,strName);
    }else{
     i = i +1;
    }
  };
  if(isFound){
    cmd_val = coboCommands.values[i];
  }else{
    cmd_val = -1;
  }
  
  return cmd_val;
}
int FindLVPSCommand(char *cmd){
  int cmd_val;
  int isFound = 0;
  int i = 0;
  int comparisonSize = 0;
  char strName[STREAM_LENGTH];
  while(!isFound && i<NLVPSCOMMANDS){
    ClearString(strName);
    GetDictionaryEntry(strName,lvpsCommands.names[i]);    
    RemoveTrailingCharacters(strName);
    comparisonSize = FindSmaller(strlen(cmd),strlen(strName));
    if(strcmp(cmd,strName)==0){
     isFound = 1;
     //printf("Here");
    }else{
     i = i +1;
    }
  };
  if(isFound){
    cmd_val = lvpsCommands.values[i];
  }else{
    cmd_val = -1;
  }
  return cmd_val;
}
int FindPDUCommand(char *cmd){
  int cmd_val;
  int isFound = 0;
  int i = 0;
  int comparisonSize = 0;
  char strName[STREAM_LENGTH];
  while(!isFound && i<NPDUCOMMANDS){
    ClearString(strName);
    GetDictionaryEntry(strName,pduCommands.names[i]);    
    RemoveTrailingCharacters(strName);
    comparisonSize = FindSmaller(strlen(cmd),strlen(strName));
    if(strcmp(cmd,strName)==0){
     isFound = 1;
     //printf("Here");
    }else{
     i = i +1;
    }
  };
  cmd_val = pduCommands.values[i];
  return cmd_val;
}
int FindSequenceCommand(char *cmd){
  int cmd_val;
  int isFound = 0;
  int i = 0;
  int comparisonSize = 0;
  char strName[STREAM_LENGTH];
  while(!isFound && i<NSEQUENCECOMMANDS){
    ClearString(strName);
    GetDictionaryEntry(strName,seqCommands.names[i]);    
    RemoveTrailingCharacters(strName);
    CleanString(strName);
    comparisonSize = FindSmaller(strlen(cmd),strlen(strName));
    //printf("%s %s\n", cmd, strName);
    fflush(stdout);
    if(strncmp(cmd,strName,comparisonSize)==0){
     isFound = 1;
     //printf("Here");
    }else{
     i = i +1;
    }
  };
  if(isFound){
    cmd_val = seqCommands.values[i];
  }else{
    cmd_val = -1;
  }
  
  return cmd_val;
}

int FindGyroCommand(char *cmd){
  int cmd_val;
  int isFound = 0;
  int i = 0;
  int comparisonSize = 0;
  char strName[STREAM_LENGTH];
  while(!isFound && i<NGYROCOMMANDS){
    ClearString(strName);
    GetDictionaryEntry(strName,gyroCommands.names[i]);    
    RemoveTrailingCharacters(strName);
    comparisonSize = FindSmaller(strlen(cmd),strlen(strName));
    if(strcmp(cmd,strName)==0){
     isFound = 1;
     //printf("Here");
    }else{
     i = i +1;
    }
  };
  if(isFound){
    cmd_val = gyroCommands.values[i];
  }else{
    cmd_val = -1;
  }

  //printf("%d\n",cmd_val);
  
  return cmd_val;
}




int SearchCommand(int id, char *cmd){
  int cmd_val;
  int isPayloadNeeded = 0;
  
  RemoveTrailingCharacters(cmd);
  switch (id){
    case QID_RC:
      cmd_val = FindRemoteControlCommand(cmd);  
      break;
    case QID_TRGB:
      cmd_val = FindTriggerBoardCommand(cmd);
      break;
    case QID_SIAB:
      cmd_val = FindSIABCommand(cmd);
      break;
    case QID_COBO:
      cmd_val = FindCoBoCommand(cmd);
      break;
    case QID_LVPS:
      cmd_val = FindLVPSCommand(cmd);
      break;
    case QID_PDU:
      cmd_val = FindPDUCommand(cmd);
      break;
    case QID_SEQ:
      cmd_val = FindSequenceCommand(cmd);
      break;
    case QID_GYRO:
      cmd_val = FindGyroCommand(cmd);
      break;
    default:
      //printf("Command not found");
      cmd_val = -1;
      break;
  };

  return cmd_val;
}

int PopulateCmdElements(char **elements, char *in){
  int whitespace = ' ';
  
  char restOfStr[STREAM_LENGTH+1] = {0};
  
  strcpy(restOfStr, in);
  RemoveTrailingCharacters(restOfStr);
  //printf("%d\n",strlen(restOfStr));
  int prevOffset = 0;
  int offset = 0;
  int i = 0;
  int isMultipleArgs = 0;
  char *placeOfWhitespace = strchr(restOfStr,whitespace);
  //printf("%s\n",placeOfWhitespace);
  //printf("%d",strlen(in));
  while(placeOfWhitespace != NULL){
    //printf("WP: %d\n",placeOfWhitespace-restOfStr);

    offset = placeOfWhitespace - restOfStr;
    //printf("\n%d",offset);
    strncpy(elements[i],restOfStr,offset);
   // printf("\n%s\n",elements[i]);
    //printf("Copy Success!");
    //fflush(stdout);
    ClearString(restOfStr);

    strncpy(restOfStr,in+prevOffset+offset+1,strlen(in)-offset-prevOffset-1);
    placeOfWhitespace = strchr(restOfStr,whitespace);
    i = i +1;
    
    //RemoveTrailingCharacters(placeOfWhitespace);
    //RemoveTrailingCharacters(restOfStr);
    //printf("%s\n",placeOfWhitespace);
    
    prevOffset = prevOffset+offset+1;
    //printf("\n%d\n",prevOffset);
  };
  //printf("\n%d\n",strlen(restOfStr));
  strncpy(elements[i],restOfStr,strlen(restOfStr));

  
  //ClearString(restOfStr);
  //ClearString(placeOfWhitespace);

  return i;
}

int CheckPayloadNeeded(int id, int cmd){
  int isPayloadNeeded = 0;
  switch (id){
    case QID_RC:
      switch (cmd){
        case UPDATE_RUNTIME:
        case GET_LOG_ID:
        case DO_TRG_SCAN:
        case DO_TRG_SCAN_SF:
        case HV_SCAN:
        case DOWNLOAD_STARLINK:
          isPayloadNeeded = 1;
          break;
        default:
          isPayloadNeeded = 0;
          break;
        };
      break;
    case QID_TRGB:
      switch(cmd){
        case CMD_PREPARE:
          isPayloadNeeded = 1;
          break;
        case CMD_CONFIG:
          isPayloadNeeded = 1;
          break;
        case CMD_READ:
        case CMD_WRITE:
          isPayloadNeeded = 1;
          break;
        default:
          isPayloadNeeded = 0;
          break;
        };
      break;
    case QID_PDU:
      switch(cmd){
        case CMD_TURN_ON:
        case CMD_TURN_OFF:
        case CMD_GET_STATE:
        case CMD_GET_USAGE:
          isPayloadNeeded = 1;
          break;
        default:
          isPayloadNeeded = 0;
          break;
      };
      break;
    case QID_COBO:
      isPayloadNeeded = 0;
      break;
    case QID_SEQ:
      switch(cmd){
        case CONFIGURE_SEQ:
          isPayloadNeeded = 1;
          break;
        default:
          isPayloadNeeded = 0;
          break;
      }
      break;
    case QID_SIAB:
      isPayloadNeeded = 1;
      break;
    case QID_GYRO:
      isPayloadNeeded = 1;
      break;
    case QID_LVPS:
      switch(cmd){
        case SET_LOOP_DELAY:
        case SET_UPDATED_HV_VALUE:
          isPayloadNeeded = 1;
          break;
        default:
          isPayloadNeeded = 0;
          break;
      };

    default:
      break;
  };

  return isPayloadNeeded;
}

unsigned short GetHighByte(int value){
  return (value&0xFF00)>>8;
}
unsigned short GetLowByte(int value){
 return value&0x00FF;
}

void ParseRemoteControlPayload(int id, int cmd, int args, char **elements, char* command){
  //char payload[RC_PAYLOAD*2] = {0};
  unsigned short start, nSteps, stepSz;
  unsigned short clientID, nLines,runtime, currentLim;
  unsigned short dwnldType;

  char timeInt = 1;

  int position;
  start = 0;
  nSteps = 0;
  stepSz = 0;
  clientID = 0;
  nLines = 0;
  runtime = 0;
  currentLim = 0;
  dwnldType = 0;



  char j = 0;
  int isFound = 0;

  char *flags[] = {"-s","-ns","-sz","-t","-pr", "-nl","-rt","-i","-ty"};
        
  
  switch(cmd){
    case UPDATE_RUNTIME:
      for(int i = 2;i<args;i+=2){
          if(strcmp(elements[i], flags[6])==0){
            runtime = atoi(elements[i+1]);
            if(runtime>0xFFFF){
              runtime=0xFFFF;
            }if(runtime<0){
              runtime = 1;
            }
        }
      }
      sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetHighByte(runtime),
                        GetLowByte(runtime),0,0,0,0,0);
      break;
    case GET_LOG_ID:
      for(int i = 2;i<args;i+=2){
          if(strcmp(elements[i], flags[4])==0){
            position = FindPositionInTextArray(elements[i+1],processID.names,NTOTALPROCESS);
            //This again looks wrong but 
            if(position > 5){
              PrintRed();
              printf("Does not have a LOG. Setting to 00");
              PrintDefault();
              position = 0;
            }
            clientID = position;
        }if(strcmp(elements[i], flags[5])==0){
          nLines = atoi(elements[i+1]);
        }
      }

        sprintf(command,"%c%c%c%c%c%c%c%c%c",
                        GetLowByte(id),GetLowByte(cmd),GetLowByte(clientID),GetHighByte(nLines),GetLowByte(nLines),0,0,0,0);
      break;
    case DOWNLOAD_STARLINK:
      if(args<3){
        printf("Invalid number of arguments for command");
      }else{
        for(int i = 2; i < args; i+=2){
          if (strcmp(elements[i],flags[8])==0){
            dwnldType=atoi(elements[i+1]);
          }
        }
      }
      sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetLowByte(dwnldType),
              0,0,0,0,0,0);
      break;
    case DO_TRG_SCAN:
    case DO_TRG_SCAN_SF:
      if(args<3*2){
        printf("Invalid number of arguments for command");
      }else{
        for(int i = 2;i<args;i+=2){
          if(strcmp(elements[i], flags[0])==0){
            start = atoi(elements[i+1]);
        }if(strcmp(elements[i], flags[1])==0){
            nSteps = atoi(elements[i+1]);
        }if(strcmp(elements[i], flags[2])==0){
            stepSz = atoi(elements[i+1]);
        }if(strcmp(elements[i], flags[3])==0){
            timeInt = atoi(elements[i+1]);
            while(!isFound && j<TIME_INTERVALS){
              if(timeInt == timeInterval.values[j]){
                timeInt = j;
                isFound = 1;
              }
              j = j+1;
            };  
          }
        }
      }
      sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetHighByte(start),
                        GetLowByte(start),GetHighByte(nSteps),GetLowByte(nSteps),GetHighByte(stepSz),
                        GetLowByte(stepSz),GetLowByte(timeInt));
      break;
    case HV_SCAN:
      for(int i = 2;i<args;i+=2){
        //printf("%s %s\n", elements[i],elements[i+1]);
        if(strcmp(elements[i], flags[0])==0){
            start = (unsigned short)(100*atof(elements[i+1]));
        }if(strcmp(elements[i], flags[1])==0){
            nSteps = atoi(elements[i+1]);
            if(nSteps>255){
              nSteps = 255;
            }
        }if(strcmp(elements[i], flags[2])==0){
            stepSz = atoi(elements[i+1])/10;
        }if(strcmp(elements[i], flags[7])==0){
            currentLim = atoi(elements[i+1]);
        }
      }
        sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetHighByte(start),
                        GetLowByte(start),GetLowByte(nSteps),GetHighByte(stepSz),
                        GetLowByte(stepSz),GetLowByte(currentLim),0);
        break;

    default:
      break;
  };

  
}

int FindPositionInNumberArray(char *needle, const int *haystack, int nMax){
  int iNeedle = atoi(needle);
  int isFound = 0;
  int counter = 0;

 while(!isFound && counter<nMax){
    if(iNeedle == haystack[counter]){
      isFound = 1;
    }else{
      counter++;
    }
  };
  if(counter>=nMax){
    counter = -1;
  }
  return counter;
}

int FindPositionInTextArray(char *needle, const char **haystack, int nMax){
  int isFound = 0;
  int counter = 0;
  int comparisonSize = 0;
  CleanString(needle);
  while(!isFound && counter < nMax){
    //comparisonSize = FindSmaller(strlen(needle),strlen(haystack[counter]));
    if(strcmp(needle,haystack[counter])==0){
      isFound = 1;
    }
    else{
      counter++;
    }
  };
  if(counter>=nMax){
    counter = -1;
  }
  return counter;
}

void ParseTriggerBoardPayload(int id, int cmd, int args, char **elements, char* command){
  //char payload[RC_PAYLOAD*2] = {0};
  unsigned short trigType,channPeriod,stretchLength,preScale;
  unsigned short hledDelay,deadtime,internalTrigRate,intLED,extLED;
  unsigned short gpsLED,enBusy, transitBusy;

  unsigned short addr, nWords;
  unsigned int wordToWrite;
  char dash = '-';
  trigType = 0;
  channPeriod = 0;
  stretchLength = 0;
  preScale = 0;
  hledDelay = 0;
  deadtime = 0;
  internalTrigRate = 0;
  intLED = 0;
  extLED = 0;
  gpsLED = 0;
  enBusy = 0;

  unsigned short byte7,byte8;

  char timeInt;

  int position;

  int j = 0;
  int isFound = 0;
  int isBad = 0;

  //printf("Arguments: %d\n", args);
  char *flags[] = {"-tt","-cp","-sl","-ps", "-hd", "-dt","-if","-il","-el","-gl","-bs","-tb","-r","-n","-w"};
      
  switch(cmd){
    case CMD_PREPARE:
     // if(args<3){
     //   printf("Invalid number of arguments for command");
     // }else{
      //printf("Hello");
      for(int i = 2; i<args; i+=2){
        if(strcmp(elements[i], flags[11])==0){
          transitBusy = atoi(elements[i+1]);
        }
      }
      sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetHighByte(transitBusy),GetLowByte(transitBusy),0,0,0,0,0);
      //    }
    break;
    case CMD_CONFIG:;
      for(int i = 2;i<args;i+=2){
        //printf("%s %d\n", elements[i], i);
        //fflush(stdout);
        if(strcmp(elements[i], flags[0])==0){
          while(strncmp(elements[i+1],&dash,1)!=0){
            position = FindPositionInTextArray(elements[i+1],triggerBoardTT.names, NTTYPES);
            if(position >= 0 ){
              trigType += triggerBoardTT.values[position];
            }
            else{
              PrintRed();
              printf("Trigger %s Type not found\n",elements[i+1]);
              PrintDefault();
            }
            i++;
          };
          i--;
        }
      
        if(strcmp(elements[i], flags[1])==0){
              position = FindPositionInNumberArray(elements[i+1],timeInterval.values, TIME_INTERVALS);
              if(position >= 0 ){
                channPeriod = position;
              }
              else{
                PrintRed();
                printf("Time Interval %s not found\n",elements[i+1]);
                PrintDefault();
                isBad = 1;
              }
            }
  
        if(strcmp(elements[i], flags[2])==0){
            position = FindPositionInNumberArray(elements[i+1],triggerBoardDC.values, NDISCCODE);
            if(position >= 0 ){
              stretchLength = position;
            }else{
              PrintRed();
              printf("Stretch Length Code %s not found\n",elements[i+1]);
              PrintDefault();
              isBad = 1;
            }
          }

        if(strcmp(elements[i], flags[3])==0){
            preScale = atoi(elements[i+1]);
        }
        if(strcmp(elements[i], flags[4])==0){
            hledDelay = atoi(elements[i+1])/10;
        }
        if(strcmp(elements[i], flags[5])==0){
           position = FindPositionInNumberArray(elements[i+1],triggerBoardDTC.values, NDEADCODE);
            if(position >= 0 ){
              deadtime = position;
            }else{
              PrintRed();
              printf("DeadTime Code %s not found\n",elements[i+1]);
              PrintDefault();
              isBad = 1;
            }
        }
        if(strcmp(elements[i], flags[6])==0){
           position = FindPositionInNumberArray(elements[i+1],triggerBoardIC.values,NTRIGRATECODE);
            if(position >= 0 ){
              internalTrigRate = position;
            }else{
              PrintRed();
              printf("Trigger Rate %s not found\n",elements[i+1]);
              PrintDefault();
              isBad = 1;
            }
        }
        if(strcmp(elements[i], flags[7])==0){
           intLED = atoi(elements[i+1]);
        }
        if(strcmp(elements[i], flags[8])==0){
           extLED = atoi(elements[i+1]);
           //printf("Ext\n");
        }
        if(strcmp(elements[i], flags[9])==0){
           gpsLED = atoi(elements[i+1]);
           //printf("GPS!\n");
        }
        if(strcmp(elements[i], flags[10])==0){
           position = FindPositionInTextArray(elements[i+1],triggerBoardBusy.names,NBUSYCODE);
            if(position >= 0 ){
              enBusy = (unsigned short) triggerBoardBusy.values[position];
            }else{
              PrintRed();
              printf("Busy value %s not found\n",elements[i+1]);
              PrintDefault();
              isBad = 1;
            }
        }
      //printf ("Flag %d Checked!\n", i);
      }
      if(trigType == 0){
        trigType = 4;
      }
      if(isBad){
        PrintRed();
        printf("Aborting command!");
        PrintDefault();
        id = 0;
        cmd = 0;
        trigType = 0;
        channPeriod = 0;
        stretchLength = 0;
        preScale = 0;
        hledDelay = 0;
        byte7 = 0;
        byte8 = 0;
      }else{
        byte7 = deadtime<<4 | internalTrigRate;
        byte8 = intLED<<7|intLED<<6|extLED<<5|extLED<<4|gpsLED<<3|gpsLED<<2|enBusy;
      
      }
      //printf ("Assigning command");
      //printf("Command: %02x%02x%02x%02x%02x%02x%02x%02x%02x",id,cmd,
      //       trigType,channPeriod,stretchLength,preScale,
      //      hledDelay,byte7,byte8);
      //fflush(stdout);
      sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),
            GetLowByte(trigType),GetLowByte(channPeriod),GetLowByte(stretchLength),
            GetLowByte(preScale),GetLowByte(hledDelay),GetLowByte(byte7),GetLowByte(byte8));
      break;
    case CMD_READ:
      for(int i = 2; i<args; i+=2){
        if(strcmp(elements[i], flags[12])==0){
          addr = atoi(elements[i+1]);
        }if(strcmp(elements[i], flags[13])==0){
          nWords = atoi(elements[i+1]);
        }
      }
      sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetHighByte(addr),GetLowByte(addr),nWords,0,0,0,0);
      break;
    case CMD_WRITE:
      for(int i = 2; i<args; i+=2){
        if(strcmp(elements[i], flags[12])==0){
          addr = atoi(elements[i+1]);
        }if(strcmp(elements[i], flags[14])==0){
          wordToWrite = atoi(elements[i+1]);
        }
      }
      sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetHighByte(addr),GetLowByte(addr),GetHighByte(wordToWrite>>16),GetLowByte(wordToWrite>>16),
          GetHighByte(wordToWrite),GetLowByte(wordToWrite),0);
      
      break;

    default:
      break;
  };

  

}

int GetSiabID(char **elements, int args){
  char *flag = "-id";

  char dash = '-';

  int siabIDs = 0;
  int idVal = 1;
  
  char *placeOfDash;
  char *temp;

  char *range [2] = {0,0};

  range[0] = (char*)malloc(ARG_LENGTH);
  ClearString(range[0]);

  range[1] = (char*)malloc(ARG_LENGTH);
  ClearString(range[1]);

  int offset;

  for(int i = 2; i < args; i+=2){
    if(strcmp(elements[i], flag) == 0){
      placeOfDash = strchr(elements[i+1],dash);

      temp = elements[i+1];

      if(placeOfDash != NULL){
        offset = placeOfDash - elements[i+1];

        //printf("%d\n%s\n",offset,temp);
        
        strncpy(range[0],temp,offset);
       //printf("%s\n%s\n",range[0],range[1]);
       // fflush(stdout);
        strncpy(range[1],temp+offset+1, strlen(temp)-offset-1);
        

      }
      else{
        range[0] = temp;
        range[1] = temp;
      }
    }
  }

  //printf("Place %d:\n", offset);
  for(int i = atoi(range[0]); i < atoi(range[1])+1; i++){
    siabIDs = idVal<<i|siabIDs;
  }

  free(range[0]);
  free(range[1]);
  return siabIDs;
}
void ParseSIABPayload (int id, int cmd, int args, char **elements, char* command){
  unsigned short byte3, byte2,byte1;
  int thres;
  int addr;

  int musicID,regAddress;

  int siabID;
 //printf("Hello\n");
 // fflush(stdout);
  siabID = GetSiabID(elements,args);
 // printf("Hello\n");
 // fflush(stdout);
  musicID = 0;
  regAddress = 0;

  addr = 0;

  thres = 0;
  byte3 = 0;
  byte2 = 0;
  byte1 = 0;
  char *flags[] = {"-m","-r","-th","-v","-w","-c","-rt"};
  //printf("%d\n",cmd);
  switch(cmd){
    case CMD_LOAD_MUSIC_CONFIG_EPRM:
    case CMD_PRINT_LOADED_MUSIC_CONFIG:
    case CMD_WRITE_LOADED_CONFIG_TO_MUSIC:
    case CMD_UPDATE_MUSIC_CONFIG:
      for (int i = 2; i <args; i++){
         if(strcmp(elements[i], flags[0])==0){
            byte3 = FindPositionInTextArray(elements[i+1],musicIDs.names,NMUSICOPTS)+1;
            printf("Bottom\n");
        }
      }
      if(byte3 == 3){
        printf("This command does not accept 'both' as parameter.\nSetting to Bottom\n");
        byte3 = 1;
      }
      break;
    case CMD_SET_TRG_THRESHOLD:;
      //char *flags[] = {"-m","-th"};
      for (int i = 2; i <args; i++){
         if(strcmp(elements[i], flags[0])==0){
            byte3 = FindPositionInTextArray(elements[i+1],musicIDs.names,NMUSICOPTS)+1;
        }
        if(strcmp(elements[i], flags[2])==0){
            thres = atoi(elements[i+1]);
            byte2 = GetHighByte(thres);
            byte1 = GetLowByte(thres);
        }
      }
      break;
    case CMD_READ_EPRM_BYTE:
    case CMD_READ_EPRM_WORD:;
      //char *flags[] = {"-r"};
      for(int i = 2; i<args; i++){
        if(strcmp(elements[i], flags[1])==0){
          addr = atoi(elements[i+1]);
          byte3 = GetHighByte(addr);
          byte2 = GetLowByte(addr);
        }  
      }
      break;
    case CMD_WRITE_EPRM_BYTE:;
      //char *flags[] = {"-r","-w"};
      for (int i = 2; i <args; i++){
         if(strcmp(elements[i], flags[1])==0){
            addr = atoi(elements[i+1]);
            byte3 = GetHighByte(addr);
            byte2 = GetLowByte(addr);
        }
        if(strcmp(elements[i], flags[4])==0){
            byte1 = (unsigned short) atoi(elements[i+1]);
        }
      }
      break;
    case CMD_WRITE_EPRM_WORD:;
      //char *flags[] = {"-r","-w"};
      for (int i = 2; i <args; i++){
        if(strcmp(elements[i], flags[1])==0){
          byte3 = atoi(elements[i+1]);
        }
        if(strcmp(elements[i], flags[4])==0){
            thres = atoi(elements[i+1]);
            byte2 = GetHighByte(addr);
            byte1 = GetLowByte(addr);
        }
      }
      break;
    case CMD_SET_MAX_CURRENT:
      //char *flags[] = {"-c"};
      for(int i = 2; i<args; i++){
        if(strcmp(elements[i], flags[5])==0){
          addr = atoi(elements[i+1]);
          byte3 = GetHighByte(addr);
          byte2 = GetLowByte(addr);
        }  
      }
      break;
    case CMD_READ_FROM_MUSIC_REGISTER:;
      //char *flags[] = {"-m","-r"};
      for(int i = 2; i<args; i++){
        if(strcmp(elements[i], flags[0])==0){
          byte3 = FindPositionInTextArray(elements[i+1],musicIDs.names,NMUSICOPTS)+1;
        }if(strcmp(elements[i], flags[1])==0){
          byte2 = atoi(elements[i+1]);
        }
        if(byte3 == 3){
          PrintCyan();
          printf("This command does not accept 'both' as parameter.\nSetting to Bottom\n");
          PrintDefault();
          byte3 = 1;
        }
      }
      break;
    case CMD_WRITE_TO_MUSIC_REGISTER:;
      //char *flags[] = {'-m','-r','-v'};
      for(int i = 2; i<args; i++){
        if(strcmp(elements[i], flags[0])==0){
          musicID = FindPositionInTextArray(elements[i+1],musicIDs.names,NMUSICOPTS)+1;
          if(musicID == 3){
            PrintCyan();
            printf("This command does not accept 'both' as parameter.\nSetting to Bottom\n");
            PrintDefault();
            musicID = 1;
        }
        }if(strcmp(elements[i], flags[1])==0){
          regAddress = atoi(elements[i+1]);
        }
        if(strcmp(elements[i], flags[3])==0){
          thres = atoi(elements[i+1]);
        }
      }
      byte3 = musicID<<6|regAddress;
      byte2 = GetHighByte(thres);
      byte1 = GetLowByte(thres);
      break;
    case CMD_READ_SINGLE_CURRENT:
    case CMD_READ_ADC_OFFSET:;
      //char *flags[] = {"-c"};
      for(int i = 2; i<args; i++){
        if(strcmp(elements[i], flags[5])==0){
          addr = atoi(elements[i+1]);
          byte3 = addr & 0x000F;
        } 
      }
      break;
    case CMD_SET_ADC_RATE:;
      //char *flags[] = {"-rt"};
      for(int i = 2; i<args; i++){
        if(strcmp(elements[i], flags[6])==0){
          thres = atoi(elements[i+1]);
          if(thres>6){
            byte3 = 6;
          }else if(thres<1){
            byte3 = 1;
          }else{
            byte3 = thres;
          }
        } 
      }
      break;
    default:
      break;
  }

  //printf("%d\n",byte3);

  sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetHighByte(siabID>>16),GetLowByte(siabID>>16),
          GetHighByte(siabID),GetLowByte(siabID),GetLowByte(cmd),GetLowByte(byte3),GetLowByte(byte2),GetLowByte(byte1));

}

void ParseGyroPayload(int id, int cmd, int args, char **elements, char* command){
  //This might look weird. But actual process id for Gyro is 0x20
  //Using 0x21 just to be able to process the command differently
  //than siab commands. At this point change to actual ID

  id = 0x20;
  //printf("Parsing %02x",GetLowByte(cmd));

  sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),0,0,0,0,GetLowByte(cmd),0,0,0);


}

void ParsePDUPayload(int id, int cmd, int args, char **elements, char* command){
  unsigned short channel;
  char *flags [] = {"-d"};
  int position;
  for(int i = 2; i<args;i+=2){
    if(strcmp(elements[i], flags[0])==0){
      position = FindPositionInTextArray(elements[i+1],pduDevices.names,NPDUDEVICES);
      
      if(position == -1){
        PrintRed();
        printf("Device not found.\n");
        PrintDefault();
        channel = 0;
        id = 0;
        cmd = 0;
      }else{
        channel = pduDevices.values[position];
      }
      //printf("%d",channel);
    }
  }

    sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetLowByte(channel),0,0,0,0,0,0);

}

void ParseSequencePayload(int id, int cmd, int args, char **elements, char* command){
  char *flags []= {"-cs"};
  unsigned short seqID = 1;
  int position;
  int comparisonSize;

  //printf("%s %s", elements[2], flags[0]);
  //fflush(stdout);
  for(int i = 2; i<args;i+=2){
    comparisonSize = FindSmaller(strlen(elements[i]),strlen(flags[0]));
    if(strcmp(elements[i], flags[0])==0){
      position = FindPositionInTextArray(elements[i+1],configSequences.names,NCONFIGSEQUENCES);
      seqID = configSequences.values[position];
    }
  }

  //printf("%d",seqID);
  //fflush(stdout);

  sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetLowByte(seqID),0,0,0,0,0,0);
}

void ParseLVPSPayload(int id, int cmd, int args, char **elements, char* command){
  char *flags []= {"-p","-hv"};
  unsigned short frequency = 0;
  float voltage = 42.0;
  unsigned short voltageInt;
  unsigned short byte4,byte3,byte2,byte1;

  switch(cmd){
    case SET_LOOP_DELAY:
      for(int i = 2; i<args;i+=2){
        if(strcmp(elements[i], flags[0])==0){
          frequency = atoi(elements[i+1]);
        }
      }
      sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetHighByte(frequency),GetLowByte(frequency),0,0,0,0,0);
      break;
    case SET_UPDATED_HV_VALUE:
      for(int i = 2; i<args;i+=2){
        if(strcmp(elements[i], flags[1])==0){
          //printf("Hl");
          voltage = atof(elements[i+1]);
        }
      }
      voltageInt = voltage*100;
     // printf("%lg\n",voltage);
      byte4 = voltageInt/1000;
      byte3 = (voltageInt%1000)/100;
      byte2 = ((voltageInt%1000)%100)/10;
      byte1 = ((voltageInt%1000)%100)%10;

      //printf("%d%d%d%d",byte4,byte3,byte2,byte1);
      sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),GetLowByte(byte4),GetLowByte(byte3),GetLowByte(byte2),GetLowByte(byte1),0,0,0);
      break;
    default:
      break;
  }
  
}

void GetPayload(int id, int cmd, int args, char **elements, char* command){
  
  switch (id){
    case QID_RC:
      ParseRemoteControlPayload(id, cmd, args, elements, command);
      break;
    case QID_TRGB:
      //printf("Parsing\n");
      ParseTriggerBoardPayload(id,cmd,args,elements,command);
      break;
    case QID_SIAB:
      ParseSIABPayload(id,cmd,args,elements,command);
      break;
    case QID_GYRO:
      ParseGyroPayload(id,cmd,args,elements,command);
      break;
    case QID_PDU:
      ParsePDUPayload(id,cmd,args,elements,command);
      break;
    case QID_SEQ:
      ParseSequencePayload(id,cmd,args,elements,command);
      break;
    case QID_LVPS:
      ParseLVPSPayload(id,cmd,args,elements,command);
      break;
    default:
      break;
  };
}

int GetUserAccess(int id,int cmd,int userType){
  int userAccess = 0;
  if(userType == 0 || userType == 1){
    userAccess = 1;
  }else{
    switch(id){
      case QID_SEQ:
        userAccess = 1;
        break;
      case QID_RC:
        switch(cmd){
          case DO_TRG_SCAN:
          case DO_TRG_SCAN_SF:
          case HV_SCAN:
          case Archive_LOGS_ID:
          case STATE_MSG_DIS_ID:
          case STATE_MSG_EN_ID:
          case DOWNLOAD_STARLINK:
            userAccess = 1;
            break;
          default:
            userAccess = 0;
            break;
        }
      case QID_SIAB:
        switch(cmd){
          case CMD_TURN_ON_HV:
          case CMD_TURN_OFF_HV:
          case CMD_SET_TRG_THRESHOLD:
            userAccess = 1;
            break;
          default:
            userAccess = 0;
            break;
        }
      case QID_LVPS:
        switch(cmd){
          case SET_UPDATED_HV_VALUE:
            userAccess = 1;
            break;
          default:
            break;
        }
      default:
        userAccess = 0;
        break;
    }
  }

  return userAccess;
}

void LogCommand(char *in, char *command,char *usr){
  time_t t = time(NULL);
  char *filename[BUFSIZE];
  char *fileHistory[BUFSIZE];

  char *log_str = (char*)malloc(BUFSIZE*5);
  char *hist_str=(char*)malloc(BUFSIZE*5);
  char *envVar = "HOME";
  struct tm tm_l = *localtime(&t);

  //printf("Here");
  //fflush(stdout);
  sprintf(filename,"%s/Desktop/CTM-%d-%02d-%02d.log",getenv(envVar),tm_l.tm_year + 1900, 
        tm_l.tm_mon +1, tm_l.tm_mday);
  sprintf(fileHistory,"%s/Desktop/.CTM.hist",getenv(envVar));
  
  FILE *fpLog = fopen(filename,"a+");
  FILE *fpHist = fopen(fileHistory,"a+");

  //printf("%s\n",filename);
  fflush(stdout);
  if(fpLog !=NULL){
    snprintf(log_str,BUFSIZE-1,"[%02d:%02d:%02d]: %s@CTM: ",
                  tm_l.tm_hour,tm_l.tm_min,tm_l.tm_sec,usr);
    //printf("%s\n",log_str);
    //fflush(stdout);
    fputs(log_str,fpLog);
    //printf("Here\n");
    //fflush(stdout);
    ClearString(log_str);
    //printf("Here\n");
    //fflush(stdout);
    sprintf(log_str,"%s\n", in);
    fputs(log_str,fpLog);
    ClearString(log_str);

    sprintf(log_str,"[%02d:%02d:%02d]: %s@CTM->GCC: ",
                  tm_l.tm_hour,tm_l.tm_min,tm_l.tm_sec,usr);
    fputs(log_str,fpLog);
    ClearString(log_str);
    for(int i = 0; i<9;i++){
            //strcommandcheck[i] = command[i];
      sprintf(log_str,"%02x",GetLowByte(command[i]));
    }
    sprintf(log_str,"%s\n",command);
    fputs(log_str,fpLog);
  }
  fclose(fpLog);

  if(fpHist != NULL){
    sprintf(hist_str,"%s\n",in);
    fputs(hist_str,fpHist);
  }
  fclose(fpHist);
  free(log_str);
  free(hist_str);
}

void SendCommand(char *in, int userType, char* usr){

  char *cmd_elements[CMDLINE_ARGS];
  char procID[20];
  char command[CMD_LENGTH_BYTES*8];
  int isPayloadNeeded = 0;
  int args = 0;
  //char *curl = (char*)malloc(STREAM_LENGTH*4);

  char errbf [CURL_ERROR_SIZE];



  //CURL settings
  CURL *curl;
  CURLcode res;

  //Message Queue
  const char* qName = "/CStoRC";
  mqd_t mq = mq_open(qName, O_WRONLY);
  if(mq == -1){
    perror("mq_open");
    //break;
  }

  size_t message_len;

  /*struct curl_slist *headers = NULL;
  headers = curl_slist_append(headers, "Content-Type:application/octet-stream");

  curl = curl_easy_init();

//  curl_easy_setopt(curl,CURLOPT_URL,"https://192.168.2.111/ct/command");
  curl_easy_setopt(curl,CURLOPT_URL,"https://122.56.27.23/ct/command");

  curl_easy_setopt(curl,CURLOPT_PORT,443L);
  curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
  curl_easy_setopt(curl, CURLOPT_USERPWD, "ct:tm4Science241!");
  curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0);
  curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0);
  curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, 9);
  curl_easy_setopt(curl, CURLOPT_VERBOSE,0);
  curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, errbf);
  //curl_easy_setopt(curl, CURLOPT_COPYPOSTFIELDS, 1);*/



  //Input Processing

  for(int i = 0; i<CMDLINE_ARGS; i++){
    cmd_elements[i] = (char*)malloc(ARG_LENGTH);
    ClearString(cmd_elements[i]);
  }
  RemoveTrailingCharacters(in);
  //printf("%s",usr);
  args = PopulateCmdElements(cmd_elements, in);

  for (int i = 0; i<args+1; i++){
     CleanOutputString(cmd_elements[i]);
  }
  // //printf("\n%s\n",procID);
  // //fflush(stdout);

  int id, cmd;
  id = SearchProcessID(cmd_elements[0]);
  if(args >= 1){
    cmd = SearchCommand(id, cmd_elements[1]);
  }else{
    //id = -1;
    cmd = -1;
    PrintGray();
    printf("Insufficent arguments\n");
    PrintDefault();
  }

  char *strcommandcheck = (char*) malloc(18);
  char *strcommand = (char*) malloc(18);
  for(int i = 0; i<9; i++){
   strcommandcheck[i] = 0;
  }
  //printf("%d",id);
  if(id != -1){
    if(cmd != -1){
      if(GetUserAccess(id,cmd,userType)){
       isPayloadNeeded = CheckPayloadNeeded(id, cmd);
       if(isPayloadNeeded){
         //printf("Payload Needed %02x\n",cmd);
         //fflush(stdout);
         GetPayload(id,cmd,args,cmd_elements,command);
       }
       else{
          //printf("Here");
          fflush(stdout);
          sprintf(command,"%c%c%c%c%c%c%c%c%c",GetLowByte(id),GetLowByte(cmd),0,0,0,0,0,0,0);
        }

        if((int)command[0] == 0){
          PrintRed();
          printf("Aborting Command\n");
          PrintDefault();
        }else{
          LogCommand(in, command,usr);
          PrintGray();
          printf("Sending Command: %s\n",in);
          printf("Message being sent to the System:");
          PrintYellow();
          for(int i = 0; i<9;i++){
            strcommandcheck[i] = command[i];
            printf("%02x",GetLowByte(strcommandcheck[i]));
          }
          sprintf(strcommand,"%02x%02x%02x%02x%02x%02x%02x%02x%02x",GetLowByte(strcommandcheck[0]),GetLowByte(strcommandcheck[1]),GetLowByte(strcommandcheck[2]),GetLowByte(strcommandcheck[3]),GetLowByte(strcommandcheck[4]),GetLowByte(strcommandcheck[5]),GetLowByte(strcommandcheck[6]),GetLowByte(strcommandcheck[7]),GetLowByte(strcommandcheck[8]));
          printf("\n");
          PrintDefault();
       
        message_len = strlen(strcommand);

        if(mq_send(mq,strcommand,message_len,0)==-1){
          perror("mq_send");
          //break;
        }else{
          PrintGray();
          printf("Command successfully sent!\n");
          PrintDefault();
        }
        if (mq_close(mq) == -1) {
          perror("mq_close");
          //return 1;
        }
       /* curl_easy_setopt(curl, CURLOPT_POSTFIELDS, command);
        res = curl_easy_perform(curl);*/

        for(int i = 0; i<CMDLINE_ARGS; i++){
          free(cmd_elements[i]);
        }

        /*if(res != CURLE_OK){
          PrintRed();
          printf("%s\n", errbf);
        }else{
          PrintGray();
          printf("Command successfully sent!\n");
          PrintDefault();
        }*/
        PrintDefault();
        }
     }
      else{
       printf("%s\nUser has no access to that command\n",KGRY);
      }
    }
    else{
      printf("%s\nCommand not found\n",KGRY);
      fflush(stdout);    
    }
  }
  else{
    printf("%s\nInvalid Subsystem\n",KGRY);
  }
  printf("%s",KDEF);
}
