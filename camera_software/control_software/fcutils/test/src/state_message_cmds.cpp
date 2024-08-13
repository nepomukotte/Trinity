#include "state_message_cmds.h"

std::string Byte_to_Str(uint32_t a)
{
    std::string status(1, static_cast<unsigned char>(a));
    return status;
}

std::string Short_to_Str(uint32_t a)
{
    std::string b1 = Byte_to_Str(((a >> 8) & 0xFF));
    std::string b2 = Byte_to_Str((a & 0xFF));
    std::string result = b1 + b2;
    return result;
}

std::string Word_to_Str(uint32_t a)
{
    std::string b1 = Byte_to_Str(((a >> 24) & 0xFF));
    std::string b2 = Byte_to_Str(((a >> 16) & 0xFF));
    std::string b3 = Byte_to_Str(((a >> 8) & 0xFF));
    std::string b4 = Byte_to_Str((a & 0xFF));
    std::string result = b1 + b2 + b3 + b4;
    return result;
}

int Parse_emon(const std::string s)
{
    std::string EMON = s.substr(1,6);
    int last = EMON.size() -1;
    while (last >= 0 && EMON[last] == ' ')
        --last;
    int emon_value = stoi(EMON.substr(0, last+1));
    return emon_value;
}

void Get_EMON(int &first_emon, int &second_emon)
{
    system("candump can0 -a > emon_tmp.txt &");
    sleep(5);
    system("sudo echo $(/bin/ps -fu $USER | awk '/candump/ && !/awk/ {print $2}' | xargs kill)");

    const std::string filename = "emon_tmp.txt";
    std::ifstream file;
    file.open(filename.c_str(), std::fstream::in);
    std::string firstline;
    if(getline(file, firstline))
    {
        int length = firstline.size();
        std::string tmp1 = firstline.substr(length-10,length);

        if((tmp1 == "'EMON   1'") || (tmp1 == "'EMON   2'"))
        {
            getline(file, firstline);
            length = firstline.size();
            tmp1 = firstline.substr(length-10,length);
            std::string secondline;
            for(int i=0; i<6; i++)
            {
                getline(file, secondline); // we skip 5 lines to get to other EMON data
            }
            std::string tmp2 = secondline.substr(length-10,length);
            first_emon = Parse_emon(tmp2);
            second_emon = Parse_emon(tmp1);
            std::cout << "1st Emon: " << first_emon << " pW" << std::endl;
            std::cout << "2nd Emon: " << second_emon << " pW" << std::endl;
        }
        else{
            std::string secondline;
            getline(file, secondline);
            std::string tmp2 = secondline.substr(length-10,length);

            first_emon = Parse_emon(tmp2);
            second_emon = Parse_emon(tmp1);
            std::cout << "1st Emon: " << first_emon << " pW" << std::endl;
            std::cout << "2nd Emon: " << second_emon << " pW" << std::endl;
        }
    }else{
        std::cout << "EMON File is empty." << std::endl;
        first_emon = 0;
        second_emon = 0;
    }
    file.close();

    system("sudo rm emon_tmp.txt");
}

uint32_t *Get_LVPS_data()
{
    static uint32_t state_msg_arr[64] = {0};
    std::string get_lvps_sample = "tail -n 2 "+LOG_DIR+"Client_LVPS.log > lvps_tmp.txt &";
    system(get_lvps_sample.c_str());
    //system("cat lvps_tmp.txt &"); // For Debugging
    usleep(100000);
    std::string filename = "lvps_tmp.txt";
    std::ifstream file;
    file.open(filename.c_str(), std::fstream::in);
    std::string lvps_sample;
    getline(file, lvps_sample);
    file.close();
    system("sudo rm lvps_tmp.txt");
    if(lvps_sample.size() != 641)
    {
        std::cout << "LVPS sample size: " << lvps_sample.size() << std::endl;
        std::cout << "Size of 1st HK Sample is wrong. Trying one more time ..." << std::endl;
        system(get_lvps_sample.c_str());
        usleep(100000);
        file.open(filename.c_str(), std::fstream::in);
        getline(file, lvps_sample);
        file.close();
        system("sudo rm lvps_tmp.txt");
        if(lvps_sample.size() != 641)
        {
            std::cout << "Size of 2nd HK Sample is wrong. Skipping this time." << std::endl;
            return state_msg_arr;
        }
    }
    for(int i=0; i<64; i++)
    {
        std::string tmp = lvps_sample.substr(i*10, 10);
        std::stringstream ss;
        ss << std::hex << tmp;
        ss >> state_msg_arr[i];
        //std::cout << state_msg_arr[i] << std::endl;   // For Debugging
    }
    return state_msg_arr;
}

std::string Get_PDU_State()
{
    system("python3 /home/cherenkov/Programs/control_software/pdu/pdu.py can0 ch9 usage > pdu_tmp.txt &");
    sleep(1);
    system("python3 /home/cherenkov/Programs/control_software/pdu/pdu.py can0 ch8 usage >> pdu_tmp.txt &");
    sleep(1);

    std::string filename = "pdu_tmp.txt";
    std::string pdu_state_str;
    std::ifstream file;
    file.open(filename.c_str(), std::fstream::in);
    std::string line1, line2, line3, line4;
    getline(file, line1);
    getline(file, line2);
    getline(file, line3);
    getline(file, line4);
    file.close();
    system("sudo rm pdu_tmp.txt");

    pdu_state_str = Short_to_Str(stoi(line1)) + Short_to_Str(stoi(line2)) + Short_to_Str(stoi(line3)) + Short_to_Str(stoi(line4));
    return pdu_state_str;
}