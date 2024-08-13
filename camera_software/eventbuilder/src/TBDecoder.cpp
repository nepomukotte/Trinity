#include <fstream>
#include <iostream>
#include "TBDecoder.h"
using namespace std;

class TBDecoder {
    private:    
    //INSTANCE VARIABLES
    EventData* pEvents; //pointer to EventData array
    int N; //number of events in file

    //--------------------------------------------------------------------------
    public:    
    //CONSTRUCTORS
    TBDecoder(int N, const char *fileName) { 
        //set number of events
        this->N = N;
        //allocate memmory and return pointer
        pEvents = new EventData[N];
        //check succesful mem allocation
        //if (pEvents == NULL) exit;
        
        //read input file
        ifstream infile;
        infile.open(fileName, ios::binary | ios::in);
        //check file opened successful
        if(!infile.is_open()) {
            cout << "Cant open file" << endl;
            infile.close();
        } else {
            //read data into allocated memmory
            infile.read((char*)pEvents, N*sizeof(EventData));
            infile.close();
        }
        unshuffleData(); //reorganizes data from 3412 to 1234, function found bellow.
    }
    
    //--------------------------------------------------------------------------
    public:
    // GETTERS - getEventNumber, getClockCount, getTriggerSource, getActiveSDB
    // Assuming event indx starts from 0;

    /**
     * Getter function for event number.
     * @param indx event index
     * @return event number
     */
    uint32_t getEventNumber(uint32_t indx){
        return pEvents[indx].eventNumber;
    }

    /**
     * Getter function for clock count.
     * @param indx event index
     * @return clock count
     */
    uint64_t getClockCount(uint32_t indx){
    	Break4Bytes tempUnion;
		tempUnion.fullWord = getWord3(indx);
		uint64_t temp64 = tempUnion.quartWord[msbClockByte];
		temp64 = (temp64<<fullWordSizeBits);
		temp64 |= pEvents[indx].word2;

        //cout<<temp64<<endl;
		return temp64;
    }

    /**
	 * Provides trigger source which is the last 6 bits of word3.
	 * @param indx The event index.
	 * @return The source of the trigger. Each trigger source corresponds to an
	 * activated bit. Returning as an unsigned int we get the following return values.
	 * If more than one trigger source is activated then each bit is simply activated.
	 * For example, bifocal and GPS => 17
	 * 1: Bifocal
	 * 2: Discriminator Test
	 * 4: Internal
	 * 8: External
	 * 16: GPS
	 * 32: LED
	 */
	uint8_t getTriggerSource(uint32_t indx){
		Break4Bytes temp;
		temp.fullWord = getWord3(indx);
		uint8_t trigger = temp.quartWord[triggerSourceByte];
		trigger = (trigger >> 2);
		return trigger;
	}

	/**
	 * Function which returns the least significant set bit of the stretched
	 * discriminator bits (words 4/5). Employs a binary search algorithm.
	 * 	Ex: sdb = 01100000 00000000 00000000 00000000 00000000 00000000 00000000 00000000
	 * 	Returns: 61
	 *
	 * Note: The behavior of this function when the stretched discriminator bits contain
	 * unexpected values is not defined (and the result meaningless). We are anticipating
	 * exactly 2 adjacent bits on the same row to be activated.
	 * 	Ex: sdb = 00000000 00000000 00000000 00000000 00000000 00000001 10000000 00000000
	 * 	Calling this function with the above SDB will return 15. This may seem correct, however,
	 * 	it is not physically accurate because the two active bits are on different rows in the
	 * 	physical photomultiplier array.
	 *
	 * @param The event index. (eventNumber - 1)
	 * @return The bit index, -1 if bifocal trigger bit is not set, -1 if no active SDBs
	 */
	uint32_t getActiveSdb(uint32_t index){
		//check if bifocal first
		/*if(!(getTriggerSource(index)&bifocalTrigger)){
			return -1;
		}

		//local variables
		uint64_t sdbBits = getSdb(index); //store stretched disciminator bits in local variable
		uint16_t activeRowVal = 0; //will store potion of SDBs with set bits
		uint8_t bitNum = 0; //records bit position

		//if there are no active bits return -1;
		if(!sdbBits){
            return -1;
		}

		//Begin Binary Search

		//Reducing from 64 bits to 32 bits
		if(sdbBits >= 0x100000000){ //any of upper 4 bits set
			sdbBits = (sdbBits >> halfNumSdbBits);
			bitNum += halfNumSdbBits;
		}

		//Reducing from 32 bits to 16 bits
		if(sdbBits >= 0x10000){ //any of upper 4 bits set
			sdbBits = (sdbBits >> quartNumSdbBits);
			bitNum += quartNumSdbBits;
		}

		//Recursively continue binary search algorithm with 16 bit backing data type.
		//Reduce from 16 -> 8 -> 4 -> 2 -> 1 bit(s).
		bitNum += recursiveHelper(eithNumSdbBits, sdbBits);

		return bitNum;*/
        // Check to see if the event is bifocal
    if(!(getTriggerSource(index)&bifocalTrigger||getTriggerSource(index)&discTestTrigger)){
        //std::cout << "This event is not bifocal." << std::endl;
        return -1;
    }

    //local variables
    bool found_pair = false;
    uint32_t bitNum = 0;                    //records the triggered bit position. starts from zero.
    uint64_t sdbBits = getSdb(index);   //store stretched disciminator bits in local variable

    //if there are no active bits return -1;
    if(!sdbBits){
        std::cout << "Disc. bits are all zero. Is the event file bad?" << std::endl;
        return -1;
    }

    // Find the bit#i, and check to see if it is non-zero and no bifocal pair has been found yet.
    // If bit#i is non-zero, then check bit#i+1, and if that is non-zero as well, then the pair is found. Retrun.
    // If the bit#i+1 is zero, then keep bit#i as temporary triggered bit and continue the for loop.
    // If by the end of for loop, no pair has been found, then temporary bit will be returned. There is no other scenario here!
    for(int i=0; i<64; i++)
    {
        if(((sdbBits >> i) & 0x01) && (!found_pair))
        {
            //std::cout<<"Music: "<<i<<" triggered"<<endl;
            if((sdbBits >> i+1) & 0x01)
            {
                //std::cout << "Found a bifocal pair at Music IDs: " << i << " and " << i+1 << std::endl;
                bitNum = i;
                found_pair = true;
            }else{
                //std::cout << "Disc. bit#" << i << " was 1, but bit#" << i+1 << " was zero!" << std::endl;
                //std::cout << "Keeping bit#" << i << " as the temporary triggered bit." << std::endl;
                bitNum = i;
            }
        }
    }

    /*if(found_pair){
        std::cout << "Found a pair. Returning bit#" << bitNum << " as the final triggered bit." << std::endl;
    }else{
        std::cout << "Did not find a pair. Returning bit#" << bitNum << " as the final triggered bit." << std::endl;
    }*/

    return bitNum;

	}

	//--------------------------------------------------------------------------
    private:
    //GETTERS
    uint32_t getWord2(uint32_t indx){
		return pEvents[indx].word2;
	}
    uint32_t getWord3(uint32_t indx){
		return pEvents[indx].word3;
	}
    uint64_t getSdb(uint32_t indx){
		return combineSdb(indx);
	}
    uint32_t getSdbLow(uint32_t indx){
        return pEvents[indx].sdbLow;
    }
    uint32_t getSdbHigh(uint32_t indx){
        return pEvents[indx].sdbHigh;
    }
    
    //--------------------------------------------------------------------------
    private:
    //SETTERS - used for unshuffling data
    void setEventNumber(uint32_t indx, uint32_t data){
        pEvents[indx].eventNumber = data;
    }
    void setWord2(uint32_t indx, uint32_t data){
        pEvents[indx].word2 = data;
    }
    void setWord3(uint32_t indx, uint32_t data){
        pEvents[indx].word3 = data;
    }
    void setSdbHigh(uint32_t indx, uint32_t data){
        pEvents[indx].sdbHigh = data;
    }
    void setSdbLow(uint32_t indx, uint32_t data){
        pEvents[indx].sdbLow = data;
    }
    
    //--------------------------------------------------------------------------
    public:
    //ADDITIONAL FUNCTIONS
    /**
    *   Frees memory allocated for the data from the current file.
    *   You should no longer use this objects pointer to the EventData array.
    *   The data is no longer reserved and therefore its behavior undefined.
    *   @param none
    *   @return none
    **/
    void removeData() {
        delete[] pEvents;
    }
    
    /**
    *   Redundancy check which confirms if index number of Event matches the
    *   event number.
    *   @param indx The event index.
    *   @return Boolean true if indices match, false if not.
    **/
    bool confirmEventNumber(uint32_t indx){
        if(indx+1 == getEventNumber(indx)){
        	return true;
        }
        return false;
    }

    //--------------------------------------------------------------------------
    private:
    // HELPER FUNCTIONS
    
    /**
	 * Helper Function for our binary search algorithm. Used once data can be
	 * stored in a 16 bit variable (once we have narrowed the search to 16 bits)
	 * @param searchSize The size of our next search space in bits.
	 * @param activeVal The current 16 bit number being searched.
	 * @return
	 */
	uint8_t recursiveHelper(uint8_t searchSize, uint16_t activeVal){
		uint8_t bitNum = 0;
		//base case
		if(searchSize == 1){
			//determine which bit is set
			if(!(activeVal & bit1)){
				bitNum -=1;
			}
			return bitNum;
		}

		//breaking N bits into N/2 bits
		if(activeVal >= (1<<searchSize)){ //any of upper 4 bits set
			activeVal = (activeVal >> searchSize);
			bitNum += searchSize;
		}
		return bitNum + recursiveHelper(searchSize/2, activeVal);
	}

    /**
    *   Initially, the bytes for the data stored in the EventData array are in
    *   the order 3421.
    *   This function calls a helper function to convert that data into the
    *   correct format of 1234.
    *   @return void
    **/
    void unshuffleData() {
    	uint32_t temp;
        for(int k=0; k<N; k++) {
        	// simply call helper function for first 3 words
            // None of those functions has any effect on the
            // stored value
           // cout<<"Event No.: "<<getEventNumber(k) << endl;
            //setEventNumber(k, shiftWord(getEventNumber(k)));
            //cout<<"Event No.: "<<getEventNumber(k) << endl;

            //cout<<"Word 2: "<<getWord2(k) << endl;
            //setWord2(k, shiftWord(getWord2(k)));
            //cout<<"Word 2: "<<getWord2(k) << endl;
            //cout<<"Word 3: "<<getWord3(k) << endl;*/
            //setWord3(k, shiftWord(getWord3(k)));
            //cout<<"Word 3: "<<getWord3(k) << endl;/*
            // call helper function and swap high and low sdb words.
            
            temp = getSdbHigh(k);
            //cout<<"High: " << getSdbHigh(k) <<endl;
            setSdbHigh(k, shiftWord(getSdbLow(k)));
            //cout<<"High: " << getSdbHigh(k) <<endl;
            //cout<<"Low: " << getSdbLow(k) <<endl;
            setSdbLow(k, shiftWord(temp));
            //cout<<"Low: " << getSdbLow(k) <<endl;
        }
    }
    
    /**
    *   Return single unsigned long for sdb. 
    *   Combines sdbHigh and sdbLow by bitshifting sdbHigh and computing bitwise or.
    *   @param index index of EventData array
    *   @return combined 64bit stretched discriminator bits.
    **/
    uint64_t combineSdb(uint32_t index){
        uint64_t temp = pEvents[index].sdbHigh;
        temp = temp << fullWordSizeBits;
        temp |= pEvents[index].sdbLow;
        
        return temp;
    }
    
    /**
    *   Converts 32 bit data stored as bytes 4,3,2,1 to correct value read as 1,2,3,4
    *   Uses union Break4Bytes to enable individual byte access.
    *   @param data 32 bit data stored as 4,3,2,1
    *   @return 32 bit data shifted to 1234
    **/
    uint32_t shiftWord(uint32_t data){
        Break4Bytes temp;
        temp.fullWord = data;
        uint32_t output;
        output = temp.quartWord[3];
        output = (output << byteSize*sizeof(char)) | temp.quartWord[2];
        output = (output << byteSize*sizeof(char)) | temp.quartWord[1];
        output = (output << byteSize*sizeof(char)) | temp.quartWord[0];
        return output;
    }
};
