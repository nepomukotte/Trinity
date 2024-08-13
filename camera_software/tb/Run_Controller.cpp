#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fcntl.h>
#include <sstream>
#include <unistd.h>
#include <thread>
#include <sys/types.h>
#include <sys/stat.h>
#include <mqueue.h>

using namespace std;

#define SERVER_QUEUE_NAME		"/server"
#define CLIENT_QUEUE_NAME_TB	"/client-trigger"
#define QUEUE_PERMISSIONS 0660
#define MAX_MESSAGES 5
#define MAX_MSG_SIZE 256
#define MSG_BUFFER_SIZE MAX_MSG_SIZE + 10

int main()
{
	std::string user_input;
	mqd_t qd_server, qd_client;
	struct mq_attr attr;
	char in_buffer [MSG_BUFFER_SIZE];

	attr.mq_flags = 0;
	attr.mq_maxmsg = MAX_MESSAGES;
	attr.mq_msgsize = MAX_MSG_SIZE;
	attr.mq_curmsgs = 0;

	if ((qd_server = mq_open (SERVER_QUEUE_NAME, O_RDONLY | O_CREAT, QUEUE_PERMISSIONS, &attr)) == -1)
	{
		perror ("Server: mq_open (server)");
		exit (1);
	}

	std::cout << "Enter the message to be sent to the Trigger Board" <<std::endl;
	while(std::getline(std::cin,user_input))
	{
		if ((qd_client = mq_open (CLIENT_QUEUE_NAME_TB, O_WRONLY)) == -1)
		{
			perror ("Server: Not able to open client queue");
			exit (1);
		}
		
		int n = user_input.length();
		char out_buffer[n+1];
		strcpy(out_buffer, user_input.c_str());

		if (mq_send (qd_client, out_buffer, strlen (out_buffer) + 1, 0) == -1)
		{
			perror ("Server: Not able to send message to client");
			continue;
		}

		cout << "Server: message sent to client.\n" << std::endl;
/*
        if (mq_receive (qd_server, in_buffer, MSG_BUFFER_SIZE, NULL) == -1) {
            perror ("Client: mq_receive");
            exit (1);
        }
*/        printf ("Server: message received from client: %s\n", in_buffer);
	}

	return 0;
}