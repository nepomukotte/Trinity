# CT Event Builer

## Pre-Requsites

Software Requirements: GET Software, ROOT and ExACT

To install GET software, please see the instruction on GET repository of SPB2 Gitlab at the link below to download and install source code.

https://gitlab.com/jem-euso/euso-spb2/ct/get


## Installation

The installation of the EventBuilder is done with Make. In the main directory /eventbuilder/

_Note: The makefile has a variable that may need to be modified. GET_DIR should point to the installation directory of the GET software._

```bash
make
```

Will run the _Makefile_ included in the repository and compile all necessary files to run the EventBuilder.

The executable _EventBuilder_ will be created and placed in the _./eventbuilder/_ directory.

## Usage

### Running
To run _EventBuilder_:

```bash
./EventBuilder [path/to/]CoboFormats.xcfg [path/to/]InputFileForAsAd0 [path/to]InputFileForAsAd1 [path/to/outputDirectory/]
```

This provides the compiled code with the necessary files to be able to run.

_CoboFormats.xcfg_ is the file describing the multi-frame format and raw data file that the CoBo generates.

_[path/to/]InputFileForAsAd0_ and _[path/to/]InputFileForAsAd1_ are the raw data files saved in the computer. The _EventBuilder_ will sort through the channels and assign them to the correct pixel according to _PixelMap.csv_. It merges both files from the AsAd into one with the name of the first input file and extension _.root_

_[path/to/outputDirectory/]_ directory were the merged file will be saved. 


