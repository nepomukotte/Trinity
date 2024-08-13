#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "anyio.h"
#include "lbp16.h"
#include "eth_boards.h"
#include "TB_Map.h"
#include "TB_commands.h"
#include <mqueue.h>
#include <signal.h>
#include <pthread.h>

#define TB_Q_NAME "/writetoTB"
#define RC_Q_NAME "/writetoRC"
#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)

bool new_cmd = false;
bool run_started = false;
char response[5] = "3030";

uint32_t* HextoInt(char* msg, int msg_size)
{
	static uint32_t cmd_array[8];
	for (int i = 0; i < msg_size; i=i+2)
	{
		char tmp[3];
		memcpy(tmp, &msg[i],2);
		tmp[3]='\0';
		cmd_array[i/2] = (int) strtol(tmp, NULL, 16);
	}
	int a = 0;	// Do not remove this line. It will cause hanging.
	return cmd_array;
}

void Parse_Run_CMD(char* msg, char *response)
{
	uint32_t* cmd_array;
	uint32_t msg_size = strlen(msg);

	if (msg_size == 16){
		cmd_array = HextoInt(msg, msg_size);
		access_board = Process_CMD(cmd_array, response, access_board);
	}else{
		printf("msg size: %d\n", msg_size);
		printf("Message Length Mismatch: Please enter correct message size.\n");
		strcpy(response, "MMMM");
		cmd_array = HextoInt(msg, 16);
		access_board = Process_CMD(cmd_array, response, access_board);
	}
}

static void read_mq(union sigval sv)
{
	char *arr;
	struct mq_attr attr;
	struct sigevent sev;
	mqd_t mq_tb = *((mqd_t *) sv.sival_ptr);

	if(mq_getattr(mq_tb, &attr) < 0)
		handle_error("mq_getattr");

	// Reregister for new messages on Q
	sev.sigev_notify = SIGEV_THREAD;
	sev.sigev_notify_function = read_mq;
	sev.sigev_notify_attributes = NULL;
	sev.sigev_value.sival_ptr = sv.sival_ptr;

	if (mq_notify(mq_tb, &sev) < 0)
		handle_error("Error during Re-register in mq_notify");

	// Read new message on the Q
	arr = malloc(attr.mq_msgsize);
	if (arr == NULL)
		handle_error("malloc");
	if(mq_receive(mq_tb, arr, attr.mq_msgsize, 0) < 0){
		handle_error("Error in mq_receive");
	}else{
		printf("Message Received:%s\n", arr);
		Parse_Run_CMD(arr, response);
	}
	new_cmd = true;
}

int main()
{
	struct sigevent sev;
	mqd_t mq_cs = mq_open(RC_Q_NAME, O_RDWR);
	if (mq_cs == (mqd_t) -1)
		handle_error("Error on CS message queue creation");

	mqd_t mq_tb = mq_open(TB_Q_NAME, O_RDONLY | O_NONBLOCK | O_CREAT, 0666, 0);
	if(mq_tb == (mqd_t) -1)
		handle_error("Error on TB message queue creation");

	// The process is registered for notification for new message on the Queue
	sev.sigev_notify = SIGEV_THREAD;
	sev.sigev_notify_function = read_mq;
	sev.sigev_notify_attributes = NULL;
	sev.sigev_value.sival_ptr = &mq_tb;

	if (mq_notify(mq_tb, &sev) < 0){
		handle_error("Error on message queue notify");
	}else{
	printf("Notify for message queue reception %s\n", TB_Q_NAME);
	}

	// emptying the queue to receive new messages
	ssize_t n = 0;
	struct mq_attr attr;
	if(mq_getattr(mq_tb, &attr) < 0)
		handle_error("Error in mq_getattr");
	char* arr = malloc(attr.mq_msgsize);
	if (arr == NULL)
		handle_error("malloc");
	while((n = mq_receive(mq_tb, arr, attr.mq_msgsize, 0) >= 0)){
		printf("Empty the Q. Msg rcvd %s\n", arr);
	}

	while(true){
		if(new_cmd){
			new_cmd = false;
			//int temp = mq_send(mq_cs, response, strlen(response), 0);
			printf("response:%s\n",response);
		}
		//memset(response, '\0', strlen(response));
		strcpy(response, "3030");
		usleep(100000);
	}

	mq_close(mq_tb);
	return 0;
}