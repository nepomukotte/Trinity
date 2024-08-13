TARGET = linux

CC = gcc
RM = rm -f
MATHLIB = -lm
OPT = -O0

#DEBUG = -g -pedantic -Wall -Wextra
#DEBUG = -g -Wall -Wextra
DEBUG = -g

ifeq ($(TARGET),linux)
	INCLUDE = -I/usr/include
	BIN = client_tb
	LIBS = -lpci -lrt
endif

CFLAGS = $(INCLUDE) $(OPT) $(DEBUG) $(MATHLIB)

objects = common.o bitfile.o hostmot2.o eeprom.o anyio.o eth_boards.o eeprom_local.o eeprom_remote.o main.o TB_commands.o

headers = eth_boards.h anyio.h hostmot2.h lbp16.h common.h eeprom.h eeprom_local.h eeprom_remote.h bitfile.h hostmot2_def.h TB_commands.h

$(BIN) all : $(objects)
	$(CC) -o $(BIN) $(objects) $(MATHLIB) $(LIBS)

main.o : main.c $(headers)
	$(CC) $(CFLAGS) -c main.c

anyio.o : anyio.c $(headers)
	$(CC) $(CFLAGS) -c anyio.c

TB_commands.o : TB_commands.c $(headers)
	$(CC) $(CFLAGS) -c TB_commands.c

eth_boards.o : eth_boards.c $(headers)
	$(CC) $(CFLAGS) -c eth_boards.c

eeprom_local.o : eeprom_local.c $(headers)
	$(CC) $(CFLAGS) -c eeprom_local.c

eeprom_remote.o : eeprom_remote.c $(headers)
	$(CC) $(CFLAGS) -c eeprom_remote.c

hostmot2.o : hostmot2.c $(headers)
	$(CC) $(CFLAGS) -c hostmot2.c

eeprom.o : eeprom.c $(headers)
	$(CC) $(CFLAGS) -c eeprom.c

bitfile.o : bitfile.c $(headers)
	$(CC) $(CFLAGS) -c bitfile.c

common.o : common.c $(headers)
	$(CC) $(CFLAGS) -c common.c

clean :
	$(RM) $(BIN) *.o

.PHONY: install
install: $(BIN)
	install --mode=0755 --owner root --group root --dir $(DESTDIR)/bin
	install --mode=0755 --owner root --group root $(BIN) $(DESTDIR)/bin