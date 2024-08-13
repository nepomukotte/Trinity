#ifndef _STATE_MSG_H
#define _STATE_MSG_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <thread>
#include <functional>
#include "constants.h"

using namespace std;

std::string Byte_to_Str(uint32_t a);
std::string Short_to_Str(uint32_t a);
std::string Word_to_Str(uint32_t a);
int Parse_emon(const std::string s);
void Get_EMON(int &first_emon, int &second_emon);
uint32_t *Get_LVPS_data();
std::string Get_PDU_State();

#endif