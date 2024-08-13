//Control Definitions
#ifndef CONTROLS_H
#define CONTROLS_H
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define STREAM_LENGTH 150
#define KEY_BACKSPACE 127
#define KEY_ENTER 10
#define KEY_EXIT 4
#define KEY_NAV 27
#define KEY_LEFT 68
#define KEY_RIGHT 67
#define KEY_UP 65
#define KEY_DOWN 66
#define KEY_TAB 9

#define KDEF  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYLW  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define KGRY  "\033[1;30m"


void Backspace();
void PrintPrompt();
void SendCommand(char *in, int userType,char *usr);
int SearchCommand(int id, char *cmd);
int SearchProcessID(char *pID);
void GetDictionaryEntry(char *dictEntry, const char* dict);
void RemoveTrailingCharacters(char* str);
void ClearString(char* str);
int PopulateCmdElements(char **elements, char *in);
void ParseRemoteControlPayload(int id, int cmd, int args, char **elements, char* command);
void GetPayload(int id, int cmd, int args, char **elements, char* command);
int CheckPayloadNeeded(int id, int cmd);

#endif