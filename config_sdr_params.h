/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#ifndef __CONFIG_SDR_PARAMS__
#define __CONFIG_SDR_PARAMS__

#include <string>
#include <cstdint>
#include <fstream>

#include <complex>
#include "logger.h"

/*******************************************************
* NAMESPACES
********************************************************/
using namespace std;
using namespace logger;


namespace config
{

/*******************************************************
* Constants
********************************************************/
static const int32_t numTotalSatellites    = 32;
static const int32_t numTotalAcqAlgorithms = 5;

static const int32_t cohrProcessTimeMs   = 2;
static const int32_t incohrProcessTimeMs = 10;

static constexpr int32_t s_acqSatellitePrnList[] = {
    /* MATLAB: S = sprintf('%i,',1:32); S(1:end) */
    1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
    17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32
};

static const double s_acqDopplerBwKHz  = 50;
static const double s_acqDopplerResHz  = 1000;
static const double s_acqThreshold     = 4;
static const int32_t s_numAcqSatellites = 8;

static const int32_t numTotalDataFiles = 10;
static const string dataInfoFileName   = "data_file_list.txt";
static const string dataInRelPath      = "..\\..\\data\\data_in\\";
static const string dataOutRelPath     = "..\\..\\data\\data_out\\";

static const double speedOfLight = 299792458;

static const string s_acqAlgoList[] = {
    "norm_acq_parcode",
    "weak_acq_dbzp"/*,
    "norm_acq_parfreq",
    "weak_acq_dbzp2",
    Add more
    */
};



/*******************************************************
* STRUCTS
********************************************************/
typedef enum MemType_e
{
    // This memory is allocated in preprocessing and 
    // is used in processing among all frames per file.
    // To be de-allocated at the end of file processing.
    PRE_PROCESS_MEM_ALLOCATED,

    // This memory is allocated in processing and 
    // is allocated in the first frame of each file
    // in the processing stage, re-used in subsequent 
    // frames processing and de-allocated at the end of 
    // file processing.
    PROCESS_MEM_ALLOCATED,

    POST_PROCESS_MEM_ALLOCATED,

    // This mean un-allocated or already de-allocated memory.
    MEM_INVALID

} MemType_t;


/**
* Signal structures
*/
typedef struct PreProcessSignals_s
{
	int32_t          optimizeOption;
	int32_t          averFactor;
	int32_t          dopplerResHz;
	int32_t          numBlocks;
	int32_t          numSamplesPerBlock;
    complex<double>* caCodeTable;
    MemType_t        caCodeTableMemType;
    int32_t          caCodeTableBytesAllocated;
	complex<double>* dopplerCplxExp;
    MemType_t        dopplerCplxExpMemType;
    int32_t          dopplerCplxExpBytesAllocated;
	complex<double>* dopplerCplxExpDelta;
    MemType_t        dopplerCplxExpDeltaMemType;
    int32_t          dopplerCplxExpDeltaBytesAllocated;
    int32_t          shiftFactor;
    int32_t          numDopplerSamples;
} PreProcessSignals_t;

typedef struct ProcessSignals_s
{
    complex<double>* ddMap;
    MemType_t        ddMapMemType;
    int32_t          ddMapBytesAllocated;
    int32_t          optimizeOption;
    int32_t          averFactor;
    int32_t          numCodeSamples;
    int32_t          numDopplerBins;
    int32_t          dopplerResHz;
    int32_t          numSamplesPerBlock;
    int32_t          numBlocks;
    int32_t          numDopplerFftBins;
} ProcessSignals_t;

typedef struct PrnProcessingResults_s
{
    double  dopplerShiftHz;
    int32_t codeDelay;
    double  peakMetric;
    double  estIfFreqHz;
    int32_t satellitePrn;

}PrnProcessingResults_t;

typedef struct PostProcessResults_s
{
    complex<double>* accumDdm;
    MemType_t        accumDdmMemType;
    int32_t          accumDdmBytesAllocated;
    int32_t          numSamplesPerPrn;
    int32_t          numAcqSatellites;
    PrnProcessingResults_t prnResults[numTotalSatellites];
} PostProcessResults_t;


typedef struct RxDataChannelMem_s
{
	bool      isInit;
	int32_t*  rxDataPerFrame;
    MemType_t rxDataPerFrameMemType;
    int32_t   rxDataPerFrameBytesAllocated;
} RxDataChannelMem_t;



/** System parameters for 
* 
*/
typedef struct SysParams_s
{
    string   acqAlgosList[numTotalAcqAlgorithms];
    int32_t  numAcqAlgos;
	int32_t  coherentProcessingTimeMS;
	int32_t  inCoherentProcessingTimeMS;
	double   startOffset;
	double   caCodeChipRateHz;
    double   minSamplingFreqHz;
	int32_t  sampleInterpOrder;
    int32_t  numberOfChannels;
    int32_t  skipNumberOfBytes;
    bool     skipAcquisition;
    int32_t  acqSatelliteList[numTotalSatellites];
    int32_t  numCfgSatellites;
    double   acqDopplerBwKhz;
    double   acqDopplerResHz;
    double   acqThreshold;
    int32_t  numAcqSatellites;
    int32_t  msToProcess;
    double   dllDampingRatio;
    double   dllNoiseBandwidth;
    double   dllCorrelatorSpacing;
    double   pllDampingRatio;
    double   pllNoiseBandwidth;
    double   navSolPeriod;
    int32_t  elevationMask;
    bool     useTropCorr;
    int32_t  truePositionE;
    int32_t  truePositionN;
    int32_t  truePositionU;
} SysParams_t;


typedef struct DataParams_s
{
    double  intermFreqHz;
    double  samplingFreqHz;
    string  dataType;
    bool    isBasebandSignal;
    int32_t totalChannels;
    int32_t selectedChannel;
} DataParams_t;


typedef struct StateParams_s
{
    string  dataPathIn;
    string  dataOutPath;
    string  fileNames[numTotalDataFiles];
    int32_t numFilesProcessed;
    int32_t numFilesToProcess;
    int32_t currFrameNum;
    int32_t numTotalFrames;
} StateParams_t;

typedef struct sdrParams_s
{
    SysParams_t   sysParams;
    StateParams_t stateParams;
    DataParams_t  dataParamsList[numTotalDataFiles];
} SdrParams_t;

/*******************************************************
* FUNCTION DECLARATIONS
********************************************************/

/**
* @brief Configures sdrParams structure parameters.
* @param [inout] p_sdrParams
* returns void
*/
void ConfigSdrParams(SdrParams_t* p_sdrParams);

}

// Parse file information here.

#endif __CONFIG_SDR_PARAMS__