/** \file T2KEvent.h
    \brief Header of the class T2KEvent.
*/
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//
//
//	RESTSoft : Software for Rare Event Searches with TPCs
//
//      T2KEvent.h
//
//      Class containing a T2K event.
//
//      History of corrections:
//
//      Dec 2009:  First concept 
//                 F.J. Iguaz
//
//_____________________________________________________________________________

#ifndef REST_T2KEvent
#define REST_T2KEvent

#include "TObject.h"
#include "TArrayI.h"

/** \class T2KEvent
    \brief Class containing a T2K event.
*/
//_____________________________________________________________________________
class T2KEvent : public TObject
{
 public:

    Int_t   fNChannels;  ///< Number of channels.
    Int_t   fNBins;      ///< Number of bins per channel.
    Int_t   fNTBins;     ///< Total number of bins.

    /// Array with all the data.
    TArrayI fData;       // [fNTBins] 

    /// Constructor.
    T2KEvent();
    /// Destructor.
    ~T2KEvent();

    // Setters.
    /// Number of channels is set.
    void SetfNChannels(int num) {
      fNChannels = num; SetfNTBins();
    }

    /// Number of bins is set.
    void SetfNBins(int num) {
      fNBins = num; SetfNTBins();
    }

    /// Total number of bins is set.
    void SetfNTBins() {
      fNTBins = fNChannels * fNBins;
      fData.Set(fNTBins);
    }

    // Getters.
    TArrayI GetChannel(int num) {
      int i;
      TArrayI channel;

      channel.Set(fNBins);

      for(i = 0; i < fNBins; i++)
	channel.fArray[i] = fData.fArray[i + fNBins * num];

      return channel;
    }

    // Setters.
    void SetPoint(int channel, int bin, double value) {
      fData.fArray[bin + fNBins * channel] = (int) value;
    }

    /// Definition of the dictionary.
    ClassDef(T2KEvent, 1);
};
//_____________________________________________________________________________
#endif
