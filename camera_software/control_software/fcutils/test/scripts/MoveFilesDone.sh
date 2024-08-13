#!/bin/bash

FILEDIR="/home/cherenkov/Documents/TemporaryData"
TBDIR="/home/cherenkov/Programs/control_software/tb/mesaflash"
DESTDIR=$1

sudo mv ${FILEDIR}/format/*AsAd0* ${DESTDIR}
sudo mv ${FILEDIR}/format/*AsAd2* ${DESTDIR}
sudo mv ${TBDIR}/TB_data* ${DESTDIR}