#include "EventBuilder.h"
#include <fstream>

/**
 * Prints out various parameters for a specified eventID in a TBDecoder object, a pointer to which
 * is stored as an EventBuilder instance variable.
 * @param TBEvents is a TBDecoder object which stores the data for N events from the TriggerBoard.
 * @param eventID is an integer specifying the event to print.
 */

void EventBuilder::logTBevent(int eventID, TBDecoder *TBEvents)
{
    eventNumber_TB = TBEvents->getEventNumber(eventID);
    eventTime_TB = TBEvents->getClockCount(eventID);
    triggerSource_TB = TBEvents->getTriggerSource(eventID);
    FirstTrigMusic_TB = TBEvents->getActiveSdb(eventID);

    // Optional output for tracking events.
    cout << "TB: Event number: " << eventNumber_TB << endl;
    // cout << "TB: Clock Count: " << eventTime_TB << endl;
    // cout << "TB: Trigger Source: " << triggerSource_TB << endl;
    /*cout << "TB: Confirm Event #: " << TBEvents->confirmEventNumber(eventID) << endl;

    if(!TBEvents->confirmEventNumber(eventID))
    {
        char errorMessage[50];
        sprintf(errorMessage, "TB event number mismatch, event %d", eventID);
        throw std::runtime_error(errorMessage);
    }*/
}

/**
 * Given x and y coordinates of a pixel, returns the pixelID. Bottom left is x=y=0.
 * @param nx x coordinate
 * @param ny y coordinate
 * @return PixelID
 */
int EventBuilder::FindPixel(int nx, int ny)
{
    int pixelID;
    pixelID = (nx*4+(ny%4)+(ny/4)*128);
    return pixelID;
}


/**
 * Given a pixel, return the x and y coords via a pointer.
 * @param pixel
 * @param nx x coordinate
 * @param ny y coordinate
 */
void EventBuilder::FindBin(int pixelID,int *nx, int *ny)
{
    int MUSIC_ID = pixelID/8;
    int MUSIC_Channel = pixelID%8;
    int MUSIC_row = MUSIC_ID/16;
    int MUSIC_column = MUSIC_ID%16;

    *ny = 4*MUSIC_row+MUSIC_Channel%4; 
    *nx = 2*MUSIC_column+MUSIC_Channel/4; 
}

/**
 * Given the ID of the first triggered music chip, returns a vector containing the sipm ID for each pixel
 * in the region of interest (ROI).
 */
vector<int> EventBuilder::FindNeighborPixels(int FirstTrigMusic)
{
    vector<int> vNeighbor;
    int nx, ny;
    int FirstTrigPix = FirstTrigMusic * 8;
    FindBin(FirstTrigPix, &nx, &ny);
    
    for (int j=ny-1; j<ny+5 && j<16; j++)
    {
        for (int i=nx-1; i<nx+5 && i<32; i++)
        {
            if((i > -1) && (j > -1))
            {
                int tmp = FindPixel(i, j);
                if((tmp > -1) && (tmp < 512))
                {
                    vNeighbor.push_back(tmp);
                }
            }
        }
    }
    std::sort (vNeighbor.begin(), vNeighbor.end());
    return vNeighbor;
}

/**
 * Function for dealing with frames not comming in event by event. Stores indices of correct order of frames in two
 * vectors (one for each asad, stored as member vars)
 */
int EventBuilder::FrameCounter(string CoBo_filename, int asadIdx)
{
    f_raw.open(CoBo_filename, std::ios::in | std::ios::binary);
    vAsAd0.clear();
    vAsAd1.clear();

    size_t frameCount = 0;
    size_t maxFrames = 5000;
    mfm::Frame::seekFrame(f_raw, 0);

    
    cout<<CoBo_filename<<endl;

    try {
        do {
            if (frameCount >= maxFrames) 
                break;
            try {
                curFrame = mfm::Frame::read(f_raw);
            }
            catch (const std::ifstream::failure & e) {
                if (f_raw.rdstate() & std::ifstream::eofbit) {
                    LOG_WARN() << "FrameCounterEOF reached.";
                    break;
                }
                else {
                    LOG_ERROR() << "Error reading frame: " << e.what();
                }
                return false;
            }
            if ((curFrame->header().frameType() !=0x1) and (curFrame->header().frameType() != 2))
            {
                LOG_DEBUG() << "Skipping frame #" << frameCount << "...";
                continue;
            }
            else if (frameCount%50 == 1)
            {
                LOG_DEBUG() << "Adding frame #" << frameCount << " to tree...";
            }
            frameCount++;
        }
        while (f_raw.good());
    }
    catch (const std::exception & e) {
        LOG_ERROR() << e.what();
    }

    cout << "Total Number of Events from AsAd#" << asadIdx << ": " << frameCount << endl;
    f_raw.close();

    return frameCount;
}


/**
 * Extract information from the frame header
 * itemCount_CoBo=4 for frame type 1 and itemCount_CoBo=2 for frame type 2
 * eventTime_CoBo = start time of frame relative to start of run
 * eventIdx_CoBo = the event to which this frame belongs in a run, starts at 0
 * asadIdx_CoBo =  asad index from which the frame originates.
 */
void EventBuilder::logCoBoEvent(std::auto_ptr<mfm::Frame> &tempFrame){
    itemCount_CoBo   = tempFrame->itemCount();
    eventTime_CoBo = tempFrame->headerField("eventTime").value<uint64_t>();
    eventIdx_CoBo  = tempFrame->headerField("eventIdx").value<uint32_t>();
    asadIdx_CoBo    = tempFrame->headerField("asadIdx").value<uint8_t>();
    cout << "EventID: " << eventIdx_CoBo << " | AsAd: " << (short)asadIdx_CoBo << endl;
}

/*
void EventBuilder::Set_Branch_Bifocal(Event *eventBifocal){
    Bifocal_tree->Branch("Events","Event",&eventBifocal,64000,0);
    /*
    Bifocal_tree->Branch("TBTime", &Bifocal_struct->tb_time, "tb_time/l");
    Bifocal_tree->Branch("CoBoTime", &Bifocal_struct->cobo_time, "cobo_time/l");
    Bifocal_tree->Branch("UNIXTime", &Bifocal_struct->unix_time, "unix_time/l");
    Bifocal_tree->Branch("ROIMusicID", Bifocal_struct->ROIMusicID, "ROI_MusicID[2]/i");
    Bifocal_tree->Branch("NofPixels", &Bifocal_struct->NofPixels, "N_Pixels/i");
    Bifocal_tree->Branch("ROIPixelsID", Bifocal_struct->ROIPixelID, "ROI_Pixel_ID[36]/i");
    Bifocal_tree->Branch("ROIPixelsSignalValue", Bifocal_struct->ROISignalValue, "ROI_Signal_Value[36][512]/i");
    Bifocal_tree->Branch("Pedestal", Bifocal_struct->Pedestal, "Pedestal[512]/F");   
}

void EventBuilder::Set_Branch_Background(Event *eventForced){
    Background_tree->Branch("Events","Event",&eventForced,64000,0);
    /*
    Background_tree->Branch("TBTime", &Background_struct->tb_time, "tb_time/l");
    Background_tree->Branch("CoBoTime", &Background_struct->cobo_time, "cobo_time/l");
    Background_tree->Branch("UNIXTime", &Background_struct->unix_time, "unix_time/l");
    Background_tree->Branch("SignalValue", Background_struct->SignalValue, "Signal_Value[512][512]/i");
}

void EventBuilder::Set_Branch_Hled(Event *eventHled){
    Hled_tree->Branch("Events","Event",&eventHled,64000,0);
    /*Hled_tree->Branch("TBTime", &Hled_struct->tb_time, "tb_time/l");
    Hled_tree->Branch("CoBoTime", &Hled_struct->cobo_time, "cobo_time/l");
    Hled_tree->Branch("UNIXTime", &Hled_struct->unix_time, "unix_time/l");
    Hled_tree->Branch("SignalValue", Hled_struct->SignalValue, "Signal_Value[512][512]/i");    
}*/
void EventBuilder::Write_TestEvent(string CoBo_filename_0, Event *eventTest, int asadIdx, int event_index, std::vector<std::vector<Int_t>> &signalValue, long startTime){
     if(event_index == 0){
        if(asadIdx == 0){
            f_raw_0.open(CoBo_filename_0, std::ios::in | std::ios::binary);
            curFrame = mfm::Frame::read(f_raw_0);
        }else{
            f_raw_1.open(CoBo_filename_0, std::ios::in | std::ios::binary);
            curFrame = mfm::Frame::read(f_raw_1);
        }
    // Any frame after the First event in CoBo file.
    }else{
        if(asadIdx == 0){
            curFrame = mfm::Frame::read(f_raw_0);
        }else{
            curFrame = mfm::Frame::read(f_raw_1);
        }
    }

    
    // ----------------------------------------------------------------------------------------------------------------
    // Store timing data in Bifocal struct
    // UNIX TIME WILL BE CHANGED IN FUTURE
    if(asadIdx == 0)
    {
        eventTest->SetCoBoTime((ULong64_t)eventTime_CoBo);
        eventTest->SetTBTime((ULong64_t)eventTime_TB+startTime);
        //eventTest->SetTBTime(time(NULL));
        
        
        eventIdx_CoBo_asad0 = eventIdx_CoBo;
    }else { // confirm ASAD1 data corresponsds to the same event as ASAD0.
        try{
            if(eventIdx_CoBo != eventIdx_CoBo_asad0){
                throw 1; //throw exception if 
            }
        }catch(int e){
            cout << "Event indices of ASAD_0 and ASAD_1 do not match!";
        }    
    }

    // PEDESTAL ----------------------------------------------------------------------------------------
    //  Pedestals are calculated in ExACT. No need for Pedestal Calculation
    // loop through each item in frame to store pedestal values.

    // vectors for storing channel and bucket indices
    std::vector<uint32_t> chanIdx = std::vector<uint32_t>(numChips, 0u);
    std::vector<uint32_t> buckIdx = std::vector<uint32_t>(numChips, 0u);


    
    // REGION OF INTEREST -----------------------------------------------------------------------------------

    // store active music IDs
    vector<int> roiMUSICs{FirstTrigMusic_TB,FirstTrigMusic_TB};
    eventTest->SetROIMusicID(roiMUSICs);

    /*
    Bifocal_struct->ROIMusicID[0] = FirstTrigMusic_TB;
    Bifocal_struct->ROIMusicID[1] = FirstTrigMusic_TB+1;
    */

    // find the pixels we would like to record for a bifocal event
    vector<int> roiPixels = FindNeighborPixels(FirstTrigMusic_TB);
    eventTest->SetROIPixelID(roiPixels);

    // read the first item in the frame
    item = curFrame->itemAt(0);
    field = item.field("");
    agetIdxField = field.bitField("agetIdx");
    sampleValueField = field.bitField("sample");
    uint32_t agetIdx = agetIdxField.value<uint32_t>();
    uint32_t sampleValue = sampleValueField.value<uint32_t>();
    pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
    signalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
    //Hled_struct->SignalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
    chanIdx[agetIdx]++;
    
    // loop through rest of the items in the frame
    for (size_t itemId=1; itemId < (((MaxNofChannels)+16)*MaxTimeBucket); ++itemId){
        //set current item
        item = curFrame->itemAt(itemId);
        field = item.field(field); // contains 2 main values, agetID and sampleValue
        agetIdxField = field.bitField(agetIdxField);
        sampleValueField = field.bitField(sampleValueField);
        agetIdx = agetIdxField.value<uint32_t>();
        sampleValue = sampleValueField.value<uint32_t>();
        
        // increment channel index and bucket index
        if (chanIdx[agetIdx] >= MaxChannelsPerAget)
            {
                chanIdx[agetIdx] = 0u;
                buckIdx[agetIdx]++;
            }
        
        if ((chanIdx[agetIdx]!=11) && (chanIdx[agetIdx]!=22) && (chanIdx[agetIdx]!=45) && (chanIdx[agetIdx]!=56)){
            pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
            signalValue[sipmID][buckIdx[agetIdx]] = sampleValue;

            //Hled_struct->SignalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
        }
        chanIdx[agetIdx]++;
    }
}

void EventBuilder::Write_BifocalEvent(string CoBo_filename_0, Event *eventBifocal, int asadIdx, int event_index, std::vector<std::vector<Int_t>> &signalValue, long startTime){
    // ----------------------------------------------------------------------------------------------------------------
    // First event in CoBo file, load first frame.
    if(event_index == 0){
        if(asadIdx == 0){
            f_raw_0.open(CoBo_filename_0, std::ios::in | std::ios::binary);
            curFrame = mfm::Frame::read(f_raw_0);
        }else{
            f_raw_1.open(CoBo_filename_0, std::ios::in | std::ios::binary);
            curFrame = mfm::Frame::read(f_raw_1);
        }
    // Any frame after the First event in CoBo file.
    }else{
        if(asadIdx == 0){
            curFrame = mfm::Frame::read(f_raw_0);
        }else{
            curFrame = mfm::Frame::read(f_raw_1);
        }
    }

    
    // ----------------------------------------------------------------------------------------------------------------
    // Store timing data in Bifocal struct
    // UNIX TIME WILL BE CHANGED IN FUTURE
    if(asadIdx == 0)
    {
        eventBifocal->SetCoBoTime((ULong64_t)eventTime_CoBo);
        eventBifocal->SetTBTime((ULong64_t)eventTime_TB+startTime);
        //eventBifocal->SetTBTime(time(NULL));
        
        /*
        Bifocal_struct->tb_time = (ULong64_t)eventTime_TB;
        Bifocal_struct->cobo_time = (ULong64_t)eventTime_CoBo;
        Bifocal_struct->unix_time = time(NULL);*/
        eventIdx_CoBo_asad0 = eventIdx_CoBo;
    }else { // confirm ASAD1 data corresponsds to the same event as ASAD0.
        try{
            if(eventIdx_CoBo != eventIdx_CoBo_asad0){
                throw 1; //throw exception if 
            }
        }catch(int e){
            cout << "Event indices of ASAD_0 and ASAD_1 do not match!";
        }    
    }

    // PEDESTAL ----------------------------------------------------------------------------------------
    //  Pedestals are calculated in ExACT. No need for Pedestal Calculation
    // loop through each item in frame to store pedestal values.

    // vectors for storing channel and bucket indices
    std::vector<uint32_t> chanIdx = std::vector<uint32_t>(numChips, 0u);
    std::vector<uint32_t> buckIdx = std::vector<uint32_t>(numChips, 0u);


    //Read the first item in the frame
    /*
    item = curFrame->itemAt(0);
    field = item.field("");
    agetIdxField = field.bitField("agetIdx");
    sampleValueField = field.bitField("sample");
    uint32_t agetIdx = agetIdxField.value<uint32_t>();
    uint32_t sampleValue = sampleValueField.value<uint32_t>();
    //fSampleValue = (float)sampleValue; //convert sample value to float for pedestal.
    pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
    */
    
    /*Bifocal_struct->Pedestal[sipmID] = fSampleValue;
    chanIdx[agetIdx]++;

    for (size_t itemId=1; itemId < (((MaxNofChannels/2)+16)*BiPedestalTimeBucket/2); ++itemId)
    {
        //set current item
        item = curFrame->itemAt(itemId);
        field = item.field(field); // contains 2 main values, agetID and sampleValue
        agetIdxField = field.bitField(agetIdxField);
        sampleValueField = field.bitField(sampleValueField);
        agetIdx = agetIdxField.value<uint32_t>();
        sampleValue = sampleValueField.value<uint32_t>();
        // increment channel index and bucket index
        if (chanIdx[agetIdx] >= MaxChannelsPerAget)
            {
                chanIdx[agetIdx] = 0u;
                buckIdx[agetIdx]++;
            }

        // Averaging over the pedestal values for each pixel
        if ((chanIdx[agetIdx]!=11) && (chanIdx[agetIdx]!=22) && (chanIdx[agetIdx]!=45) && (chanIdx[agetIdx]!=56)){
            fSampleValue = (float)sampleValue; //convert sample value to float for pedestal.
            pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
            if(buckIdx[agetIdx] == 0){
                Bifocal_struct->Pedestal[sipmID] = fSampleValue;
            }else if(buckIdx[agetIdx] == (BiPedestalTimeBucket-1)){
                Bifocal_struct->Pedestal[sipmID] = (Bifocal_struct->Pedestal[sipmID] + fSampleValue) / BiPedestalTimeBucket;
            }else{
                Bifocal_struct->Pedestal[sipmID] = Bifocal_struct->Pedestal[sipmID] + fSampleValue;
            }
        }
        chanIdx[agetIdx]++;
    }*/
    
    // REGION OF INTEREST -----------------------------------------------------------------------------------

    // store active music IDs
    vector<int> roiMUSICs{FirstTrigMusic_TB, FirstTrigMusic_TB+1};
    eventBifocal->SetROIMusicID(roiMUSICs);

    /*
    Bifocal_struct->ROIMusicID[0] = FirstTrigMusic_TB;
    Bifocal_struct->ROIMusicID[1] = FirstTrigMusic_TB+1;
    */

    // find the pixels we would like to record for a bifocal event
    vector<int> roiPixels = FindNeighborPixels(FirstTrigMusic_TB);
    eventBifocal->SetROIPixelID(roiPixels);
    /*for(int k=0; k<roiPixels.size(); k++){
        Bifocal_struct->ROIPixelID[k] = roiPixels.at(k);
    }
    Bifocal_struct->NofPixels = roiPixels.size();*/

    // vectors for storing channel and bucket indices
    // chanIdx = std::vector<uint32_t>(numChips, 0u);
    // buckIdx = std::vector<uint32_t>(numChips, 0u);

    /*uint16_t desiredAsad[ROIMaxNofChannels];
    uint16_t desiredAget[ROIMaxNofChannels];
    uint16_t desiredChannel[ROIMaxNofChannels];*/


    //std::vector<std::vector<Int_t>> signalValue;
    //signalValue = std::vector<std::vector<Int_t>>(MaxNofChannels, std::vector<Int_t>(MaxTimeBucket));

    // read the first item in the frame
    item = curFrame->itemAt(0);
    field = item.field("");
    agetIdxField = field.bitField("agetIdx");
    sampleValueField = field.bitField("sample");
    uint32_t agetIdx = agetIdxField.value<uint32_t>();
    uint32_t sampleValue = sampleValueField.value<uint32_t>();
    pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
    signalValue[sipmID][buckIdx[agetIdx]] = sampleValue;

    //Hled_struct->SignalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
    chanIdx[agetIdx]++;
    
    // loop through rest of the items in the frame
    for (size_t itemId=1; itemId < (((MaxNofChannels)+16)*MaxTimeBucket); ++itemId){
        //set current item
        item = curFrame->itemAt(itemId);
        field = item.field(field); // contains 2 main values, agetID and sampleValue
        agetIdxField = field.bitField(agetIdxField);
        sampleValueField = field.bitField(sampleValueField);
        agetIdx = agetIdxField.value<uint32_t>();
        sampleValue = sampleValueField.value<uint32_t>();
        
        // increment channel index and bucket index
        if (chanIdx[agetIdx] >= MaxChannelsPerAget)
            {
                chanIdx[agetIdx] = 0u;
                buckIdx[agetIdx]++;
            }
        
        if ((chanIdx[agetIdx]!=11) && (chanIdx[agetIdx]!=22) && (chanIdx[agetIdx]!=45) && (chanIdx[agetIdx]!=56)){
            pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
            signalValue[sipmID][buckIdx[agetIdx]] = sampleValue;

            //Hled_struct->SignalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
        }
        chanIdx[agetIdx]++;
    }

    //eventBifocal->SetSignalValue(signalValue);


    /*          OBSOLETE           */

    /*
    Want to save traces of all pixels in memory
    ExACT will only keep the traces for ROI pixels
    when preparing file download
    */

    /*for(int k=0; k<roiPixels.size(); k++){
        pixelMapDecoder(pixelMapArray, roiPixels.at(k), &desiredAsad[k], &desiredAget[k], &desiredChannel[k]);
        // cout << roiPixels.at(k) << ", " << desiredAsad[k] << ", " << desiredAget[k] << ", " << desiredChannel[k] << endl;
    }

    // read the first item in the frame
    item = curFrame->itemAt(0);
    field = item.field("");
    agetIdxField = field.bitField("agetIdx");
    sampleValueField = field.bitField("sample");
    agetIdx = agetIdxField.value<uint32_t>();
    sampleValue = sampleValueField.value<uint32_t>();
    pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
    for(int k=0; k<roiPixels.size(); k++){
        if(desiredAsad[k] == asadIdx){
            if(desiredAget[k] == agetIdx){
                if(desiredChannel[k] == chanIdx[agetIdx]){
                    Bifocal_struct->ROISignalValue[k][buckIdx[agetIdx]] = sampleValue;
                    break;
                }
            }
        }
    }
    chanIdx[agetIdx]++;

    for (size_t itemId=1; itemId < (((MaxNofChannels/2)+16)*MaxTimeBucket/2); ++itemId){
        //set current item
        item = curFrame->itemAt(itemId);
        field = item.field(field); // contains 2 main values, agetID and sampleValue
        agetIdxField = field.bitField(agetIdxField);
        sampleValueField = field.bitField(sampleValueField);
        agetIdx = agetIdxField.value<uint32_t>();
        sampleValue = sampleValueField.value<uint32_t>();

        // increment channel index and bucket index
        if (chanIdx[agetIdx] >= MaxChannelsPerAget)
            {
                chanIdx[agetIdx] = 0u;
                buckIdx[agetIdx]++;
            }

        if ((chanIdx[agetIdx]!=11) && (chanIdx[agetIdx]!=22) && (chanIdx[agetIdx]!=45) && (chanIdx[agetIdx]!=56)){
            pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
            for(int k=0; k<roiPixels.size(); k++){
                if(desiredAsad[k] == asadIdx){
                    if(desiredAget[k] == agetIdx){
                        if(desiredChannel[k] == chanIdx[agetIdx]){
                            Bifocal_struct->ROISignalValue[k][buckIdx[agetIdx]] = sampleValue;
                            break;
                        }
                    }
                }
            }
        }
        chanIdx[agetIdx]++;
        
    }*/
}

void EventBuilder::Write_BackgroundEvent(string CoBo_filename_0, Event *eventForced, int asadIdx, int event_index, std::vector<std::vector<Int_t>> &signalValue, long startTime){
    // ----------------------------------------------------------------------------------------------------------------
    // First event in CoBo file, load first frame.
    if(event_index == 0){
        if(asadIdx == 0){
            f_raw_0.open(CoBo_filename_0, std::ios::in | std::ios::binary);
            curFrame = mfm::Frame::read(f_raw_0);
        }else{
            f_raw_1.open(CoBo_filename_0, std::ios::in | std::ios::binary);
            curFrame = mfm::Frame::read(f_raw_1);
        }
    // Any frame after the First event in CoBo file.
    }else{
        if(asadIdx == 0){
            curFrame = mfm::Frame::read(f_raw_0);
        }else{
            curFrame = mfm::Frame::read(f_raw_1);
        }
    }

    // -----------------------------------------------
    // Debugging - Print Cobo header file data
    //-----------------------------------------------
    cout << "--------------------------------" << endl;
    cout << "From Cobo Header:" << endl;
    cout << "ASAD: " << curFrame->headerField("asadIdx").value<uint32_t>() << endl;
    cout << "Event Indx: " << curFrame->headerField("eventIdx").value<uint32_t>() << endl;
    cout << "Event Time: " << curFrame->headerField("eventTime").value<uint64_t>() << endl;
    cout << "--------------------------------" << endl;

    // ----------------------------------------------------------------------------------------------------------------
    // Store timing data in background struct
    // UNIX TIME WILL BE CHANGED IN THE FUTURE
    if(asadIdx == 0){

        eventForced->SetCoBoTime((ULong64_t)eventTime_CoBo);
        eventForced->SetTBTime((ULong64_t)eventTime_TB + startTime);
        //eventForced->SetTBTime(time(NULL));
        /*
        Background_struct->tb_time = eventTime_TB;
        Background_struct->cobo_time = eventTime_CoBo;
        Background_struct->unix_time = time(NULL);
        */
        eventIdx_CoBo_asad0 = eventIdx_CoBo;
    }else { // Confirm ASAD1 data corresponds to the same event as ASAD0
        try{
            if(eventIdx_CoBo != eventIdx_CoBo_asad0){
                throw 1; //throw exception if 
            }
        }catch(int e){
            cout << "Event indices of ASAD_0 and ASAD_1 do not match!";
        }
    }
    // ----------------------------------------------------------------------------------------------------------------
    // Store item data in background struct

    // vectors for storing channel and bucket indices
    std::vector<uint32_t> chanIdx = std::vector<uint32_t>(numChips, 0u);
    std::vector<uint32_t> buckIdx = std::vector<uint32_t>(numChips, 0u);
    
    //std::vector<std::vector<Int_t>> signalValue;
    //signalValue = std::vector<std::vector<Int_t>>(MaxNofChannels, std::vector<Int_t>(MaxTimeBucket));


    // read the first item in the frame
    item = curFrame->itemAt(0);
    field = item.field("");
    agetIdxField = field.bitField("agetIdx");
    sampleValueField = field.bitField("sample");
    uint32_t agetIdx = agetIdxField.value<uint32_t>();
    uint32_t sampleValue = sampleValueField.value<uint32_t>();
    pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
    
    signalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
    //Background_struct->SignalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
    chanIdx[agetIdx]++;

    // loop through each item in frame
    for (size_t itemId=1; itemId < (((MaxNofChannels)+16)*MaxTimeBucket); ++itemId){
        //set current item
        item = curFrame->itemAt(itemId);
        field = item.field(field); // contains 2 main values, agetID and sampleValue
        agetIdxField = field.bitField(agetIdxField);
        sampleValueField = field.bitField(sampleValueField);
        agetIdx = agetIdxField.value<uint32_t>();
        sampleValue = sampleValueField.value<uint32_t>();
        
        // increment channel index and bucket index
        if (chanIdx[agetIdx] >= MaxChannelsPerAget)
            {
                chanIdx[agetIdx] = 0u;
                buckIdx[agetIdx]++;
            }
        
        if ((chanIdx[agetIdx]!=11) && (chanIdx[agetIdx]!=22) && (chanIdx[agetIdx]!=45) && (chanIdx[agetIdx]!=56)){
            pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
            signalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
            //Background_struct->SignalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
        }
        chanIdx[agetIdx]++;
    }

    //eventForced->SetSignalValue(signalValue);
}

void EventBuilder::Write_HledEvent(string CoBo_filename_0, Event *eventHled, int asadIdx, int event_index, std::vector<std::vector<Int_t>> &signalValue, long startTime){
    // ----------------------------------------------------------------------------------------------------------------
    // First event in CoBo file, load first frame.
    if(event_index == 0){
        if(asadIdx == 0){
            f_raw_0.open(CoBo_filename_0, std::ios::in | std::ios::binary);
            curFrame = mfm::Frame::read(f_raw_0);
        }else{
            f_raw_1.open(CoBo_filename_0, std::ios::in | std::ios::binary);
            curFrame = mfm::Frame::read(f_raw_1);
        }
    // Any frame after the First event in CoBo file.
    }else{
        if(asadIdx == 0){
            curFrame = mfm::Frame::read(f_raw_0);
        }else{
            curFrame = mfm::Frame::read(f_raw_1);
        }
    }
    
    // -----------------------------------------------
    // Debugging - Print Cobo header file data
    //-----------------------------------------------
    /*cout << "--------------------------------" << endl;
    cout << "From Cobo Header:" << endl;
    cout << "ASAD: " << curFrame->headerField("asadIdx").value<uint32_t>() << endl;
    cout << "Event Indx: " << curFrame->headerField("eventIdx").value<uint32_t>() << endl;
    cout << "Event Time: " << curFrame->headerField("eventTime").value<uint64_t>() << endl;
    cout << "--------------------------------" << endl;
*/
    // ------------------------------------Filename: ./R----------------------------------------------------------------------------
    // Store cobo header data and TB Time (from TB header data) in Hled struct
    // UNIX TIME WILL BE CHANGED IN THE FUTURE.

    if(asadIdx == 0){
        eventHled->SetCoBoTime((ULong64_t)eventTime_CoBo);
        eventHled->SetTBTime((ULong64_t)eventTime_TB+startTime);
        //cout<<(ULong64_t)eventTime_TB+startTime<<endl;
        //cout<<startTime<<endl;
        cout<<eventTime_TB<<endl;
        //eventHled->SetTBTime(time(NULL));
        /*
        Hled_struct->tb_time = eventTime_TB;
        Hled_struct->cobo_time = eventTime_CoBo;
        Hled_struct->unix_time = time(NULL);
        */
        eventIdx_CoBo_asad0 = eventIdx_CoBo;
    }else { // confirm ASAD1 data corresponds to the same event as ASAD0.
        try{
            if(eventIdx_CoBo != eventIdx_CoBo_asad0){
                throw 1;
            }
        }catch(int e){
            cout << "Event indices of ASAD_0 and ASAD_1 do not match!";
        }
    }

    // ----------------------------------------------------------------------------------------------------------------
    // Store item data in Hled struct
    // vectors for storing channel and bucket indices
    std::vector<uint32_t> chanIdx = std::vector<uint32_t>(numChips, 0u);
    std::vector<uint32_t> buckIdx = std::vector<uint32_t>(numChips, 0u);

    //std::vector<std::vector<Int_t>> signalValue;
    //signalValue = std::vector<std::vector<Int_t>>(MaxNofChannels, std::vector<Int_t>(MaxTimeBucket));

    // read the first item in the frame
    item = curFrame->itemAt(0);
    field = item.field("");
    agetIdxField = field.bitField("agetIdx");
    sampleValueField = field.bitField("sample");
    uint32_t agetIdx = agetIdxField.value<uint32_t>();
    uint32_t sampleValue = sampleValueField.value<uint32_t>();
    pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
    signalValue[sipmID][buckIdx[agetIdx]] = sampleValue;

    //Hled_struct->SignalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
    chanIdx[agetIdx]++;

    // loop through rest of the items in the frame
    for (size_t itemId=1; itemId < (((MaxNofChannels)+16)*MaxTimeBucket); ++itemId){
        //set current item
        item = curFrame->itemAt(itemId);
        field = item.field(field); // contains 2 main values, agetID and sampleValue
        agetIdxField = field.bitField(agetIdxField);
        sampleValueField = field.bitField(sampleValueField);
        agetIdx = agetIdxField.value<uint32_t>();
        sampleValue = sampleValueField.value<uint32_t>();
        
        // increment channel index and bucket index
        if (chanIdx[agetIdx] >= MaxChannelsPerAget)
            {
                chanIdx[agetIdx] = 0u;
                buckIdx[agetIdx]++;
            }
        
        if ((chanIdx[agetIdx]!=11) && (chanIdx[agetIdx]!=22) && (chanIdx[agetIdx]!=45) && (chanIdx[agetIdx]!=56)){
            pixelMapEncoder(pixelMapArray, &sipmID, asadIdx, agetIdx, chanIdx[agetIdx]);
           /* if(sipmID==0){
                cout<<sipmID<<" "<<buckIdx[agetIdx]<<" "<< sampleValue<<endl;
            }*/
            
            signalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
            //Hled_struct->SignalValue[sipmID][buckIdx[agetIdx]] = sampleValue;
        }

        
        chanIdx[agetIdx]++;
    }

    //eventHled->SetSignalValue(signalValue);
}
/*
void EventBuilder::Write_Tree_to_file(){
    Bifocal_tree->Write();
    Background_tree->Write();
    Hled_tree->Write();
}*/


//Event Builder Constructor
EventBuilder::EventBuilder(string formats){
    {
//	cout<<"ME"<<endl;
        mfm::FrameDictionary::instance().addFormats(formats.c_str());

        // create TTrees for storing data and create ouput tree file.
//	cout<<formats<<endl;        

        // load pixelMap csv file and pixelMap
        char* eventEnv = std::getenv("EVENTBUILDER_DIR");
        string path(eventEnv);
        string fileName = path+"/Pixel_Map.csv";
//	cout<<fileName<<endl;
        pixelMapArray = openPixelMap(fileName);
    }
}

//Event Builder Destructor
EventBuilder::~EventBuilder(){
    closePixelMap(pixelMapArray);
}

void EventBuilder::mainFlow(string filename_0, string filename_tb, string outDir)
{
    struct std::tm  t = {};

    std::string CoBo_filename_0 = filename_0 + ".graw";
    //std::string CoBo_filename_1 = filename_1 + ".graw";
    std::string Out_filename = filename_0 + ".root";
    std::string TB_filename = filename_tb + ".bin";

    std::string TB_trimmedName = TB_filename.substr(TB_filename.find_last_of("/\\")+1);
    TB_trimmedName = TB_trimmedName.substr(0,TB_trimmedName.find(".bin")-5);

    strptime(TB_trimmedName.c_str(), "TB_data_%y_%m_%d_%H_%M_%S",&t);

    long seconds = mktime(&t)*1000*1000*100; //In 10s of ns from Jan 1 1970

    Out_filename = Out_filename.substr(Out_filename.find_last_of("/\\")+1);
    Out_filename = outDir+"/"+Out_filename;

    cout<<"R"<<endl;
    /*
    STRUCT_Bifocal Bifocal_event;
    STRUCT_Background Background_event;
    STRUCT_Hled Hled_event;
    */
    /* Here we switch to Event DataType*/
    

    TFile *f = new TFile(Out_filename.c_str(), "RECREATE");

    TTree *Bifocal_tree = Bifocal_tree = new TTree("BiFocal", "Bifocal Triggers");
    TTree *Background_tree = new TTree("Forced", "Forced Triggers");;
    TTree *Hled_tree = new TTree("HLED", "HLED Triggers");
    TTree *Test_tree = new TTree("Test", "Test Triggers");

    Event *eventBifocal = new Event();
    Event *eventHled = new Event();
    Event *eventForced = new Event();
    Event *eventTest = new Event();

    Bifocal_tree->Branch("Events","Event",&eventBifocal,64000,0);
    Background_tree->Branch("Events","Event",&eventForced,64000,0);
    Hled_tree->Branch("Events","Event",&eventHled,64000,0);
    Test_tree->Branch("Events","Event",&eventTest,64000,0);


/*
    Set_Branch_Bifocal(Bifocal_tree, eventBifocal);
    Set_Branch_Hled(Hled_eventHled);
    Set_Branch_Background(eventForced);
*/
    /*Set_Branch_Bifocal(&Bifocal_event);
    Set_Branch_Background(&Background_event);
    Set_Branch_Hled(&Hled_event);*/

    // Open CoBo files and count Number of events for each AsAd
    int Total_NofEvents = 0;
    int AsAd0_NofEvents = FrameCounter(CoBo_filename_0,0);
    Total_NofEvents = AsAd0_NofEvents;
    //int AsAd2_NofEvents = FrameCounter(CoBo_filename_1,2);

    /*if(AsAd0_NofEvents == AsAd2_NofEvents) {
        Total_NofEvents = AsAd0_NofEvents;
        cout << "NofEvents on both AsAd files match." << endl;
    }
    else{
        if(AsAd0_NofEvents > AsAd2_NofEvents){
        Total_NofEvents = AsAd2_NofEvents;
        cout << "NofEvents on AsAd0 file is bigger. Going with AsAd2_NofEvents." << endl;
        }else{
            Total_NofEvents = AsAd0_NofEvents;
            cout << "NofEvents on AsAd2 file is bigger. Going with AsAd0_NofEvents." << endl;
        }
    }*/

    // create TBDecoder object to read trigger board data
    TBDecoder *TBEvents = new TBDecoder(Total_NofEvents, TB_filename.c_str());
    int hledCounter, biFCounter;
    hledCounter = 0;
    biFCounter = 0;
    std::vector<std::vector<Int_t>> signalValue;
    signalValue = std::vector<std::vector<Int_t>>(MaxNofChannels, std::vector<Int_t>(MaxTimeBucket));

    // loop through each event in the data (2 frames per event)
    for (int event_index=0; event_index < Total_NofEvents; event_index++)
    {
        cout << "-----------------------------------------" << endl;
        // set TB instance variables, including triggerSource_TB
        try{
            logTBevent(event_index, TBEvents);
        }catch(int e){
            cout << "TB event number mismatch, event index " << e << endl;
            return;
        }

        //triggerSource_TB = HLED;
        switch (triggerSource_TB){
            case(Bifocal):
                f->cd();
                //cout<<"MUSIC: "<<FirstTrigMusic_TB<<endl;
                if(FirstTrigMusic_TB<=63){
                    cout << "Event " << event_index << " | BIFOCAL TRIGGER" << endl;  
                    Write_BifocalEvent(CoBo_filename_0, eventBifocal, ASAD0, event_index, signalValue, seconds);
                    //Write_BifocalEvent(CoBo_filename_1, eventBifocal, ASAD1, event_index, signalValue, seconds);               
                    eventBifocal->SetSignalValue(signalValue);
                    Bifocal_tree->Fill(); 
                    biFCounter ++;
                }else{
                    cout << "Event " << event_index << " | CORRUPT TRIGGER" << endl;
                    //f->cd();
                    Write_BackgroundEvent(CoBo_filename_0, eventForced, ASAD0, event_index,signalValue,seconds);
                    //Write_BackgroundEvent(CoBo_filename_1, eventForced, ASAD1, event_index,signalValue,seconds);
                    //f->cd();
                    eventForced->SetSignalValue(signalValue);
                    Background_tree->Fill();
                }      
                break;
            case(Internal):
                cout << "Event " << event_index << " | INTERNAL TRIGGER" << endl;           
                f->cd();
                Write_BackgroundEvent(CoBo_filename_0, eventForced, ASAD0, event_index,signalValue,seconds);
                //Write_BackgroundEvent(CoBo_filename_1, eventForced, ASAD1, event_index,signalValue,seconds);
                
                eventForced->SetSignalValue(signalValue);
                Background_tree->Fill();
                break;
            case(DiscTest):
                f->cd();
                //if(FirstTrigMusic_TB<=63){
                    cout << "Event " << event_index << " | TEST TRIGGER" << endl;           
                    
                    Write_TestEvent(CoBo_filename_0, eventTest, ASAD0, event_index,signalValue,seconds);
                   // Write_TestEvent(CoBo_filename_1, eventForced, ASAD1, event_index,signalValue,seconds);
                    eventTest->SetSignalValue(signalValue);
                    Test_tree->Fill();
                /*}else{
                    cout << "Event " << event_index << " | CORRUPT TRIGGER" << endl;
                    //f->cd();
                    Write_BackgroundEvent(CoBo_filename_0, eventForced, ASAD0, event_index,signalValue,seconds);
                    Write_BackgroundEvent(CoBo_filename_1, eventForced, ASAD1, event_index,signalValue,seconds);
                    //f->cd();
                    eventForced->SetSignalValue(signalValue);
                    Background_tree->Fill();
                }*/
                break;
            case(HLED):
                cout << "Event " << event_index << " | HLED TRIGGER" << endl;
                f->cd();
                Write_HledEvent(CoBo_filename_0, eventHled, ASAD0, event_index,signalValue,seconds);
                //Write_HledEvent(CoBo_filename_1, eventHled, ASAD1, event_index,signalValue,seconds);  // we pass asad1 here, even if 2nd asad is connected to asad2 or asad3. this is just for calculating the pixel mapping
                //f->cd();
                eventHled->SetSignalValue(signalValue);
                Hled_tree->Fill();
                hledCounter++;
                break;
            default:
                cout << "Event " << event_index << " | CORRUPT TRIGGER" << endl;
                f->cd();
                Write_BackgroundEvent(CoBo_filename_0, eventForced, ASAD0, event_index,signalValue,seconds);
               // Write_BackgroundEvent(CoBo_filename_1, eventForced, ASAD1, event_index,signalValue,seconds);
                //f->cd();
                eventForced->SetSignalValue(signalValue);
                Background_tree->Fill();
                break;
        }
    }

    cout<<"BiFocal Events: "<< biFCounter <<" HLED Events: "<<hledCounter<<endl;
     
    // delete TBEvent object
    //delete TBEvents;

    // Close CoBo files, oppened within FrameCounter function above
    f_raw_0.close();
    //f_raw_1.close();

    // create output tree file, write all three trees to that file, close the file
    
    // Bifocal_tree->Write();
    // Background_tree->Write();
    // Hled_tree->Write();
    f->cd();
    f->Write();
    //Write_Tree_to_file();
    f->Close();
}
