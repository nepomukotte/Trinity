#ifndef __SMB_MAIN_H__
#define __SMB_MAIN_H__


/*
This file contains device specific #defines
Support for new devices can be added by adding a new entry below.
See application note documentation for details.
*/


#if defined(__ATmega48__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif defined(__ATmega8__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK

#elif defined(__ATmega88__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif defined(__ATmega16__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK

#elif defined(__ATmega163__)
#define TIMER1_OVF_VECT    TIMER1_OVF1_vect
#define TWI_VECT          TWSI_vect
#define TIMSK1            TIMSK

#elif defined(__ATmega168__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif defined(__ATmega168__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif defined(__ATmega32__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWSI_vect
#define TIMSK1            TIMSK

#elif (defined(__ATmega323__))
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWSI_vect
#define TIMSK1            TIMSK

#elif defined(__ATmega406__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif (defined(__ATmega64__))
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK

#elif defined(__ATmega640__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif defined(__ATmega128__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK

#elif (defined(__ATmega328__) || defined(__ATmega328p__) || defined(__AVR_ATmega328P__))
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
//#define TIMSK1            TIMSK1

#elif defined(__ATmega1280__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif (defined(__ATmega1281__) || defined(__ATmeg1281__)) // Misspelling in ioavr.h
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif defined(__AT90CAN128__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif defined(__ATmega2560__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#elif defined(__ATmega2561__)
#define TIMER1_OVF_VECT    TIMER1_OVF_vect
#define TWI_VECT          TWI_vect
#define TIMSK1            TIMSK1

#else
#error "Device not supported by the SMBus driver! This is either because the selected device does not have a TWI module, or support for the device has not been added to 'device_specific.h' yet. Please consult application note AVR316 for information on how to add new devices."
#endif


#endif
