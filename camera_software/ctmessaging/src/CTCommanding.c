#include <termios.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <sys/ioctl.h>
#include <time.h>
#include <curl/curl.h>
#include "CTM_ConstantValues.h"

#include "CTM_Controls.h"
struct termios  orig_termios;

void disableRawMode() {
  tcsetattr(STDIN_FILENO, TCSAFLUSH, &orig_termios);
}

void enableRawMode() {
  tcgetattr(STDIN_FILENO, &orig_termios);
  atexit(disableRawMode);
  struct termios raw = orig_termios;
  raw.c_lflag &= ~(ECHO|ICANON);
  tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);
}

void enablePassword(){
  tcgetattr(STDIN_FILENO, &orig_termios);
  struct termios raw = orig_termios;
  raw.c_lflag &= ~(ECHO);
  tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);
}

void ClearInputLine(char *in){
  memset(&in[0], 0, STREAM_LENGTH);
}

void CleanOutputString(char *in){
  int i = 0;
  while(isalnum(in[i]) || in[i]==' ' || in[i] == '-' || in[i] == '.'){
    i = i+1;
  }
  *(in+i) = '\0';

}
void PrintLines(int x){
  for(int i = 0; i<x; i++){
    printf("/");
    fflush(stdout);
  }
}

void PrintBlanks(int x){
  for(int i = 0; i<x; i++){
    printf(" ");
    fflush(stdout);
  }
}

void PrintNewLine(){
  printf("\n");
  fflush(stdout);
}

void PrintPrompt(){
  printf("CTM> ");
  fflush(stdout);
}

void PrintSquares(int x){
  for(int i = 0; i<x; i++){
    printf("\u2588");
    fflush(stdout);
  }
}
void PrintLetters(int w, int h, int pad){

  /*int minLetterW = 16;
  int scale = w/minLetterW;*/
  
  
  // Line 0 //
  PrintBlanks(pad);
  PrintBlanks(6);
  PrintSquares(4);
  PrintBlanks(6);
 
  PrintBlanks(pad);
  
  PrintBlanks(1);
  PrintSquares(14);
  PrintBlanks(1);
 
  PrintBlanks(pad);
  
  PrintSquares(4);
  PrintBlanks(8);
  PrintSquares(4);
  PrintNewLine();

  // Line 1
  PrintBlanks(pad);
  PrintBlanks(4);
  PrintSquares(8);
  PrintBlanks(4);
 
  PrintBlanks(pad);
  
  PrintBlanks(7);
  PrintSquares(2);
  PrintBlanks(7);
 
  PrintBlanks(pad);
  
  PrintSquares(5);
  PrintBlanks(6);
  PrintSquares(5);
  PrintNewLine();

  // Line 2
  PrintBlanks(pad);
  PrintBlanks(3);
  PrintSquares(3);
  PrintBlanks(5);
  PrintSquares(2);
  PrintBlanks(3);
 
  PrintBlanks(pad);
  
  PrintBlanks(7);
  PrintSquares(2);
  PrintBlanks(7);
 
  PrintBlanks(pad);
  
  PrintSquares(5);
  PrintBlanks(6);
  PrintSquares(5);
  PrintNewLine();

  //Line 3
  PrintBlanks(pad);
  PrintBlanks(2);
  PrintSquares(3);
  PrintBlanks(7);
  //PrintSquares(2);
  PrintBlanks(4);
 
  PrintBlanks(pad);
  
  PrintBlanks(7);
  PrintSquares(2);
  PrintBlanks(7);
 
  PrintBlanks(pad);
  
  PrintSquares(2);
  PrintBlanks(1);
  PrintSquares(3);
  PrintBlanks(4);
  PrintSquares(3);
  PrintBlanks(1);
  PrintSquares(2);

  PrintNewLine();

  //Line 4

  PrintBlanks(pad);
  PrintBlanks(1);
  PrintSquares(3);
  PrintBlanks(9);
  //PrintSquares(1);
  PrintBlanks(3);
 
  PrintBlanks(pad);
  
  PrintBlanks(7);
  PrintSquares(2);
  PrintBlanks(7);
 
  PrintBlanks(pad);
  
  PrintSquares(2);
  PrintBlanks(2);
  PrintSquares(2);
  PrintBlanks(4);
  PrintSquares(2);
  PrintBlanks(2);
  PrintSquares(2);
  PrintNewLine();

  //Line5
  PrintBlanks(pad);
  PrintBlanks(2);
  PrintSquares(3);
  PrintBlanks(7);
  //PrintSquares(2);
  PrintBlanks(4);
 
  PrintBlanks(pad);
  
  PrintBlanks(7);
  PrintSquares(2);
  PrintBlanks(7);
 
  PrintBlanks(pad);
  
  PrintSquares(2);
  PrintBlanks(2);
  PrintSquares(3);
  PrintBlanks(2);
  PrintSquares(3);
  PrintBlanks(2);
  PrintSquares(2);
  PrintNewLine();

  //Line 6
  PrintBlanks(pad);
  PrintBlanks(3);
  PrintSquares(3);
  PrintBlanks(5);
  PrintSquares(2);
  PrintBlanks(3);
 
  PrintBlanks(pad);
  
  PrintBlanks(7);
  PrintSquares(2);
  PrintBlanks(7);
 
  PrintBlanks(pad);
  
  PrintSquares(2);
  PrintBlanks(3);
  PrintSquares(6);
  PrintBlanks(3);
  PrintSquares(2);
  PrintNewLine();

  //Line 7
  PrintBlanks(pad);
  PrintBlanks(4);
  PrintSquares(8);
  PrintBlanks(4);
 
  PrintBlanks(pad);
  
  PrintBlanks(7);
  PrintSquares(2);
  PrintBlanks(7);
 
  PrintBlanks(pad);
  
  PrintSquares(2);
  PrintBlanks(4);
  PrintSquares(4);
  PrintBlanks(4);
  PrintSquares(2);
  PrintNewLine();

   //Line 8
  PrintBlanks(pad);
  PrintBlanks(6);
  PrintSquares(4);
  PrintBlanks(6);
 
  PrintBlanks(pad);
  
  PrintBlanks(7);
  PrintSquares(2);
  PrintBlanks(7);
 
  PrintBlanks(pad);
  
  PrintSquares(2);
  PrintBlanks(5);
  PrintSquares(2);
  PrintBlanks(5);
  PrintSquares(2);
  PrintNewLine();


}

void PrintLogo(){
  struct winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ,&w);

  int minLetterW = 16;
  int minLetterH = 20;

  int maxLetterW = 48;
  int maxLetterH = 60;

  int nLetters = 3;
  int spaceBetween = 5;
  int pad = 5;
  int letterW, letterH;

  /*if(w.ws_col>=(3*maxLetterH+2*spaceBetween+2*pad) && w.ws_row>=maxLetterH){
    letterW = maxLetterW;
    letterH = maxLetterH;
  }else{
    letterW = minLetterW;
    letterH = minLetterH;
  }*/

  letterW = minLetterW;
  letterH = minLetterH;

  PrintSquares(3*minLetterW+pad*4);
  PrintNewLine();
  PrintSquares(3*minLetterW+pad*4);
  PrintNewLine();
  PrintBlanks(3*minLetterW+pad*4);
  PrintNewLine();
  PrintLetters(letterW,letterH,pad);


  PrintBlanks(3*minLetterW+pad*4);
  PrintNewLine();
  PrintSquares(3*minLetterW+pad*4);
  PrintNewLine();
  PrintSquares(3*minLetterW+pad*4);
  PrintNewLine();

}
void PrintDashes(int x){
  for(int i = 0; i<x; i++){
    printf("-");
    fflush(stdout);
  }
}
void PrintWelcomeMessage(char *usr){
  int minLetterW = 16;
  int pad = 5;
  PrintDashes(minLetterW*3+pad*4);
  PrintNewLine();
  printf("Welcome %s%s%s to CTM\n",KBLU,usr,KDEF);
  fflush(stdout);

  PrintBlanks(minLetterW*3+pad*4);
  PrintNewLine();

  printf("Usage of the system is as follows:\n");
  printf("%sSubsystem %sCommand %s[options]",KCYN,KGRN,KRED);

  PrintNewLine();
  printf("%s",KDEF);
  PrintDashes(minLetterW*3+pad*4);
  PrintNewLine();
  
  PrintBlanks(minLetterW*3+pad*4);
  PrintNewLine();
}

int Login(char *usr){
  char u[20],pwd[20];
  char *uReg[20],*uPwd[20], *uType[10];
  int isFound = 0;
  int userType = -1;
  printf("User: ");
  scanf("%s",u);
  enablePassword();
  printf("Password: ");
  scanf("%s",pwd);
  disableRawMode();

  //printf("\n%s\n%s",u,pwd);
  //fflush(stdout);

  FILE *f = fopen("./.login/.users","r");
  if(f == NULL){
    printf("Users file not found!");
  }else{
    while(fscanf(f,"%s %s %s",uReg,uPwd,uType)==3  && !isFound){  
      if(strcmp(uReg,u)==0&&strcmp(uPwd,pwd)==0){
        isFound = 1;
        userType = atoi(uType);
        //usr = u;
      }else{
        userType = -1;
      }
    }
  }
  for(int i = 0; i<strlen(u);i++){
    usr[i] = u[i];
  }
  
  //printf("%s",usr);
  return userType;
}

int LoadHistory(){
  char *envVar = "HOME";
  char *fileHistory [BUFSIZE];
  sprintf(fileHistory,"%s/Desktop/.CTM.hist",getenv(envVar));

  FILE *fpHist = fopen(fileHistory, "r");

  int nLines = 0;
  if(fpHist != NULL){
    while(fgets(fileHistory,BUFSIZE,fpHist)!=NULL){
      nLines++;
    }
  }
  return nLines;
}

int FetchHistory(int lineID, char *in){
  char *envVar = "HOME";
  char fileHistory [BUFSIZE];
  sprintf(fileHistory,"%s/Desktop/.CTM.hist",getenv(envVar));

  FILE *fpHist = fopen(fileHistory, "r");
  int isFound = 0;
  int nLines = 0;
  if(fpHist != NULL){
    while(fgets(fileHistory,BUFSIZE,fpHist)!=NULL || !isFound){
      nLines++;
      if(nLines == lineID){
        isFound = 1;
        Clear();

        for(int i = 0;i<STREAM_LENGTH;i++){
          Backspace();
        }
        PrintPrompt();
        ClearInputLine(in);
        //printf("%d",strlen(fileHistory));
        for(int i = 0; i<strlen(fileHistory); i++){
          in[i] = fileHistory[i];
          //printf("%c",fileHistory[i]);
        }
        //printf("%s",(in));
        CleanOutputString(in);
        for(int i = 0; i<strlen(in);i++){
          printf("%c",in[i]);
          fflush(stdout);
        }
        
      }
    }
  }
  fclose(fpHist);
  return strlen(in);
}

void AutoCompleteCommand(char *in, int userType){
  char *cmd_elements[CMDLINE_ARGS];
  int args = 0;
  for(int i = 0; i<CMDLINE_ARGS; i++){
    cmd_elements[i] = (char*)malloc(ARG_LENGTH);
    ClearString(cmd_elements[i]);
  }
  RemoveTrailingCharacters(in);
  args = PopulateCmdElements(cmd_elements, in);
}

int main(int argc, char* argv) {
  int userType = -1;
  int pwdAttempts = 0;
  char usr[20];
  do{
    //userType = Login(usr);
    userType=3;
    if(userType == -1){
      printf("\nInvalid username/password combination.\n");
      fflush(stdout);
    }
    if(userType == 3){
      PrintCyan();
      printf("\nHi Operator! Please type your name: ");
      PrintDefault();
      fflush(stdout);
      scanf("%s",usr);
    }
    //printf("%d",userType);
    //printf("%s",usr)
    pwdAttempts++;
  }while(userType == -1 && pwdAttempts < 3);
  userType=0;
  
  //printf("%s",usr);
  fflush(stdout);
  PrintNewLine();
  int nLinesHist = LoadHistory();
  int lineCounter = nLinesHist+1;
  if(userType != -1){

    enableRawMode();
    
    char c;
    char *input=(char*)malloc(STREAM_LENGTH);
    int ctr = 0;
    int navFlag = 0;
    

    //char *line b= 0;

    printf("%s",KGRY);
    PrintLogo();
    PrintNewLine();
    printf("%s",KDEF);
    PrintWelcomeMessage(usr);
    PrintPrompt();

    //printf("%s",usr);

    c=getchar();
  
    while (read(STDIN_FILENO, &c, 1) == 1 && (int)c!=KEY_EXIT){
      //printf("%d",(int)c);
      if (iscntrl(c)||c==KEY_NAV||c==KEY_TAB) {
        switch (c){
          case KEY_BACKSPACE:
            if(*(input) != '\0' && ctr != 0){
              if(ctr<strlen(input)){
  
                for(int i = ctr; i<strlen(input);i++){
                  input[i-1] = input[i];
                }
                Left(strlen(input)+1);
                input[strlen(input)] = '\0';
                
                Clear();
                Right(strlen(input)-ctr);
                for(int i = 0;i<STREAM_LENGTH;i++){
                  Backspace();
                }
                PrintPrompt();
                for(int i = 0; i<strlen(input); i++){
                  printf("%c",input[i]);
                  fflush(stdout);
                }
                Backspace();
                input[strlen(input)-1] = '\0';
                Left(strlen(input)-ctr+1);            
    
                ctr = ctr - 1;
              }else{
                Backspace();
                input[ctr-1]='\0';
                ctr--;
              }
            }
            break;
          case KEY_ENTER:
            if((*input) != 0){
              printf("\n");
              RemoveTrailingCharacters(input);
              CleanOutputString(input);
              
              SendCommand(input, userType,usr);
              PrintPrompt();
              fflush(stdout);

              ClearInputLine(input);
              //printf("Cleared");
              //fflush(stdout);
              //memset(input, 0, STREAM_LENGTH);
              ctr = 0;
              nLinesHist = LoadHistory();
              lineCounter = nLinesHist+1;
//              free(input);
//              input = (char*)malloc(STREAM_LENGTH);
            }else{
              printf("\n");
              PrintPrompt();
              fflush(stdout);
              memset(input, 0, STREAM_LENGTH);
//              free(input);
 //             input = (char*)malloc(STREAM_LENGTH);
              ctr = 0;
            }
            
            break;
          case KEY_NAV:
            c=getchar();
            //printf("2%d\n",c);
            //fflush(stdout);
            c=getchar();
            //printf("%d\n",c);
            //fflush(stdout);
            //mc=getchar();
  
            switch(c){
              case KEY_UP:
                //navFlag = 1;
                //printf("Up\n");
                Right(strlen(input)-ctr);
                if(lineCounter!=0){
                  lineCounter--;
                }
                ctr = FetchHistory(lineCounter,input);
                //printf("%d\n",ctr);
                
                
                break; 
              case KEY_DOWN:
                //navFlag = 1;
                //printf("Down\n");
                if(lineCounter<nLinesHist){
                  lineCounter++;
                }
                Right(strlen(input)-ctr);
                ctr = FetchHistory(lineCounter,input);
                break; 
              case KEY_LEFT:
                if(ctr>0){
                  Left(1);
                  ctr = ctr - 1;
                }
                break;
              case KEY_RIGHT:
                if(ctr<strlen(input) && (*input != 0)){
                  Right(1);
                  ctr++;
                }
                break; 
              default:
                break; 
            }
            break;
          case KEY_TAB:
            printf("Tab\n");
            break;
                
          default:
            printf("Key combination not understood\n %d",c);
            break;
        }
      }else{
        if(ctr<strlen(input) && ctr != 0){
          for(int i = strlen(input); i>ctr;i--){
            input[i] = input[i-1];
          }
          input[ctr] = (char)c;
          for(int i = ctr; i<strlen(input); i++){
            printf("%c",input[i]);
            fflush(stdout);
          }
          ctr++;
          Left(strlen(input)-ctr);
        }else{
          input[ctr] = (char)c;
          ctr = ctr +1;
          printf("%c",c);
          fflush(stdout);
        }
      }
        
    }
    printf("\n");
    }else{
      printf("\nToo many attempts\n");
    }

  return 0;
}
