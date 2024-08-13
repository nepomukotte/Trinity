#ifndef __COMMON_H
#define __COMMON_H

#include <pci/pci.h>
#define LO_BYTE(x) ((x) & 0xFF)
#define HI_BYTE(x) (((x) & 0xFF00) >> 8)

void sleep_ns(u64 nanoseconds);

#endif