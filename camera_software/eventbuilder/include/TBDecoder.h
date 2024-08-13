/**
*   Trigger Board Event Decoder Header File
*   Version: 1.0
*   Date: 5/18/2021
*   Author: Robert Pierrard
**/

#define wordSizeBytes 4
#define doubleWordSizeBytes 8
#define halfWordSizeBits 16
#define fullWordSizeBits 32
#define byteSize 8
#define numSdbRows 4
#define msbClockByte 0
#define triggerSourceByte 3

#define bifocalTrigger 0x01
#define discTestTrigger 0x02
#define halfNumSdbBits 32
#define quartNumSdbBits 16
#define eithNumSdbBits 8
#define bit1 0x0002

/**
 * Structure for storing data of 1 event consisting of 5 32bit words.
 * -----------------------------------------
 * Usage:
 * Used by TBDecoder object to store data.
 * Memory is allocated on the heap for an array of EventData structs.
 * All data in input file is then saved at one time to the eventData array.
 * Data can then be accessed in the following format: eventArray[N].eventNumber
 * -----------------------------------------
 * Note:
 * In current application, raw data is initially stored, processed, and
 * re-stored in a more readable format using a helper function.
 */
struct EventData {
	uint32_t eventNumber;
	uint32_t word2;
	uint32_t word3;
	//sdb was purposefull split into high and low bits rather than using long
	//memory allocation was problematic with variable size of long.
	uint32_t sdbHigh;
	uint32_t sdbLow;
};

/** Union for accessing any byte of 4 byte word.
 * 	--------------------------------------------
*	How to access individual bytes:
* 	- declare union and set value using fullWord member.
* 	- access individual bytes by calling quartWord with desired index.
**/
union Break4Bytes {
	uint32_t fullWord;
	uint8_t quartWord[wordSizeBytes];
};
