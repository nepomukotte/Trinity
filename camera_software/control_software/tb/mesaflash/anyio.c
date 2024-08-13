#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "anyio.h"
#include "eeprom.h"
#include "bitfile.h"
#include "eth_boards.h"

board_t boards[MAX_BOARDS];
int boards_count;

int anyio_init(board_access_t *access) {
    int ret;

    if (access == NULL)
        return -1;
    if (access->eth == 1) {
        ret = eth_boards_init(access);
        if (ret != 0)
            return ret;
    }
    return 0;
}

void anyio_cleanup(board_access_t *access) {
    if (access == NULL)
        return;
    if (access->eth == 1)
        eth_boards_cleanup(access);
}

void anyio_scan(board_access_t *access) {
    if (access == NULL)
        return;
    if (access->eth == 1)
        eth_boards_scan(access);
}

board_t *anyio_get_dev(board_access_t *access) {
    int i;
    board_t *board = NULL;

    if (access == NULL)
        return 0;
    for (i = 0; i < boards_count; i++) {
        if (strncmp(access->device_name, boards[i].llio.board_name, strlen(access->device_name)) == 0) {
            board = &boards[i];
            return board;
        }
    }
    return NULL;
}

int anyio_list_dev(board_access_t *access) {
    int i;

    if (access == NULL)
        i = 0;

    else{
        for (i = 0; i < boards_count; i++) {
        board_t *board = &boards[i];

        board->print_info(board);
        }
    }
    
    return i;
}

int anyio_dev_write_flash(board_t *board, char *bitfile_name, int fallback_flag) {
    if (board->llio.write_flash != NULL) {
        u32 addr = board->flash_start_address;

        if (fallback_flag == 1) {
            addr = FALLBACK_ADDRESS;
        }
        board->llio.write_flash(&(board->llio), bitfile_name, addr);
    } else {
        printf("ERROR: Board %s doesn't support flash writing.\n", board->llio.board_name);
        return -EINVAL;
    }
    return 0;
}

int anyio_dev_verify_flash(board_t *board, char *bitfile_name, int fallback_flag) {
    if (board->llio.verify_flash != NULL) {
        u32 addr = board->flash_start_address;

        if (fallback_flag == 1) {
            addr = FALLBACK_ADDRESS;
        }
        board->llio.verify_flash(&(board->llio), bitfile_name, addr);
    } else {
        printf("ERROR: Board %s doesn't support flash verification.\n", board->llio.board_name);
        return -EINVAL;
    }
    return 0;
}

int anyio_dev_program_fpga(board_t *board, char *bitfile_name) {
    if (board->llio.reset != NULL) {
        board->llio.reset(&(board->llio));
    } else {
        printf("ERROR: Board %s doesn't support FPGA resetting.\n", board->llio.board_name);
        return -EINVAL;
    }
    if (board->llio.program_fpga != NULL) {
        board->llio.program_fpga(&(board->llio), bitfile_name);
    } else {
        printf("ERROR: Board %s doesn't support FPGA programming.\n", board->llio.board_name);
        return -EINVAL;
    }
    return 0;
}

int anyio_dev_send_packet(board_t *board, char *lbp16_send_packet_data) {
    u8 packet[512];
    u8 *pch = lbp16_send_packet_data;
    u32 *ptr = (u32 *) packet;
    int i, recv;

    for (i = 0; i < 512, i < strlen(lbp16_send_packet_data); i++, pch += 2) {
        char s[3] = {*pch, *(pch + 1), 0};
        packet[i] = strtol(s, NULL, 16) & 0xFF;
    }
    eth_socket_send_packet(&packet, i/2);
    recv = eth_socket_recv_packet(&packet, 512);
    for (i = 0; i < recv; i++)
        printf("%02X", packet[i]);
    printf("\n");
    return 0;
}

void anyio_bitfile_print_info(char *bitfile_name) {
    FILE *fp;
    char part_name[32];

    fp = fopen(bitfile_name, "rb");
    if (fp == NULL) {
        printf("Can't open file %s: %s\n", bitfile_name, strerror(errno));
        return;
    }
    print_bitfile_header(fp, (char*) &part_name);
    fclose(fp);
}
