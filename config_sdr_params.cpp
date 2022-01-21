/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include "config_sdr_params.h"

using namespace std;

namespace config {

    /**
 * @brief Configures sdrParams structure parameters.
 * @param [inout] p_sdrParams
 * returns void
*/
    void ConfigSdrParams(SdrParams_t* p_sdrParams)
    {

        /* System Parameters */
        SysParams_t* const p_sysParams = &(p_sdrParams->sysParams);
        StateParams_t* const p_stateParams = &(p_sdrParams->stateParams);
        DataParams_t* const p_dataParams = &(p_sdrParams->dataParamsList[0]);

        p_sysParams->numAcqAlgos = sizeof(s_acqAlgoList) / sizeof(string);
        for (int32_t algoIdx = 0; algoIdx < p_sysParams->numAcqAlgos; algoIdx++)
        {
            p_sysParams->acqAlgosList[algoIdx].assign(s_acqAlgoList[algoIdx].c_str());
        }

        p_sysParams->coherentProcessingTimeMS = cohrProcessTimeMs;
        p_sysParams->inCoherentProcessingTimeMS = incohrProcessTimeMs;
        if (cohrProcessTimeMs > 10)
        {
            p_sysParams->coherentProcessingTimeMS = 10;
        }

        if (p_sysParams->inCoherentProcessingTimeMS <
            p_sysParams->coherentProcessingTimeMS)
        {
            ASSERT_WMSG("Incoherent processing time must be "
                "greater than coherent processing time.");
        }

        p_sysParams->startOffset       = 68.025f;
        p_sysParams->caCodeChipRateHz  = 1.023e6f;
        p_sysParams->minSamplingFreqHz = p_sysParams->caCodeChipRateHz * 2;
        p_sysParams->sampleInterpOrder = 2;
        p_sysParams->numberOfChannels  = 8;
        p_sysParams->skipNumberOfBytes = 0;
        p_sysParams->skipAcquisition   = 0;

        copy(begin(s_acqSatellitePrnList), end(s_acqSatellitePrnList),
            begin(p_sysParams->acqSatelliteList));

        p_sysParams->numCfgSatellites = (int32_t)size(s_acqSatellitePrnList);
        p_sysParams->acqDopplerBwKhz  = s_acqDopplerBwKHz;
        p_sysParams->acqDopplerResHz  = s_acqDopplerResHz;
        p_sysParams->acqThreshold     = s_acqThreshold;
        p_sysParams->numAcqSatellites = s_numAcqSatellites;

        p_sysParams->msToProcess          = 36;
        p_sysParams->dllDampingRatio      = 0.7f;
        p_sysParams->dllNoiseBandwidth    = 2;
        p_sysParams->dllCorrelatorSpacing = 0.5f;
        p_sysParams->pllDampingRatio      = 0.7f;
        p_sysParams->pllNoiseBandwidth    = 25;
        p_sysParams->navSolPeriod         = 500;
        p_sysParams->elevationMask        = 10;
        p_sysParams->useTropCorr          = 1;
        p_sysParams->truePositionE        = -1;
        p_sysParams->truePositionN        = -1;
        p_sysParams->truePositionU        = -1;

        /* State Parameters */
        p_stateParams->dataPathIn.assign(dataInRelPath);
        p_stateParams->dataOutPath.assign(dataOutRelPath);

        string dataInPath = p_stateParams->dataPathIn;
        ifstream fileStream(dataInPath.append(dataInfoFileName).c_str());

        int32_t fileCount = 0;
        if (fileStream)
        {
            for (string fileName; getline(fileStream, fileName); )
            {
                p_stateParams->fileNames[fileCount++].append(fileName);
            }
        }
        else
        {
            ASSERT_WMSG("Data info file not found!");
        }

        p_stateParams->numFilesProcessed = 0;
        p_stateParams->numFilesToProcess = fileCount;

        p_stateParams->currFrameNum = 0;
        p_stateParams->numTotalFrames = p_sysParams->inCoherentProcessingTimeMS /
            p_sysParams->coherentProcessingTimeMS;


        /* Data Parameters */
        for (int32_t fileIdx = 0; fileIdx < p_stateParams->numFilesToProcess; fileIdx++)
        {
            DataParams_t* const p_dataParamsFile = &p_dataParams[fileIdx];
            if (!p_stateParams->fileNames[fileIdx].compare("Woodbine_47a"))
            {
                p_dataParamsFile->intermFreqHz     = 10e6;
                p_dataParamsFile->samplingFreqHz   = 40e6;
                p_dataParamsFile->dataType         = "int16";
                p_dataParamsFile->isBasebandSignal = 1;
                p_dataParamsFile->totalChannels    = 1;
                p_dataParamsFile->selectedChannel  = 1;
            }
            else if (!p_stateParams->fileNames[fileIdx].compare("GPS_and_GIOVE_A-NN-fs16_3676-if4_1304.bin"))
            {
                p_dataParamsFile->intermFreqHz     = 4.1304e6;
                p_dataParamsFile->samplingFreqHz   = 16.367e6;
                p_dataParamsFile->dataType         = "int8";
                p_dataParamsFile->isBasebandSignal = 0;
                p_dataParamsFile->totalChannels    = 1;
                p_dataParamsFile->selectedChannel  = 1;
            }
            else if (!p_stateParams->fileNames[fileIdx].compare("GPSdata-DiscreteComponents-fs38_192-if9_55.bin"))
            {
                p_dataParamsFile->intermFreqHz     = 9.548e6;
                p_dataParamsFile->samplingFreqHz   = 38.192e6;
                p_dataParamsFile->dataType         = "int8";
                p_dataParamsFile->isBasebandSignal = 0;
                p_dataParamsFile->totalChannels    = 1;
                p_dataParamsFile->selectedChannel  = 1;
            }
            else if (!p_stateParams->fileNames[fileIdx].compare("NTLab_Bands_GPS_GLONASS_L12.bin"))
            {
                p_dataParamsFile->intermFreqHz     = (1590000000 - 1575420000);
                p_dataParamsFile->samplingFreqHz   = 53e6f;
                p_dataParamsFile->dataType         = "ubit2";
                p_dataParamsFile->isBasebandSignal = 0;
                p_dataParamsFile->totalChannels    = 4;
                p_dataParamsFile->selectedChannel  = 4;
            }
            else
            {
                ASSERT_WMSG("Listed input data file is no supported by SDR. Add support in code!");
            }
        }

    }

}
