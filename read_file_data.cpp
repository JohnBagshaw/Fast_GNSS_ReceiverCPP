/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include "cfg/logger.h"
#include "pre_process/read_file_data.h"

namespace processing
{
    void ReadDataFile(
        SdrParams_t& sdrParams,
        RxDataChannelMem_t* p_rxDataMem,
        int32_t numRxDataChannels
    )
    {

        int32_t fileIdx = sdrParams.stateParams.numFilesProcessed;
        int32_t skipNumberOfBytes = sdrParams.sysParams.skipNumberOfBytes;

        string fileName = sdrParams.stateParams.fileNames[fileIdx];
        string dataType = sdrParams.dataParamsList[fileIdx].dataType;
        double samplingFreqHz = sdrParams.dataParamsList[fileIdx].samplingFreqHz;
        int32_t selectedChannel = sdrParams.dataParamsList[fileIdx].selectedChannel == -1 ?
                                  0 : (sdrParams.dataParamsList[fileIdx].selectedChannel - 1);
        int32_t channelSeperation = sdrParams.dataParamsList[fileIdx].totalChannels;


        // Read input data 
        string fileFullName = sdrParams.stateParams.dataPathIn.append(fileName);
        FILE* fid = fopen(fileFullName.c_str(), "rb");
        if (fid == NULL)
        {
            ASSERT_WMSG("Input data file could not be read.");
        }
        else
        {

            // File opened successfully, read the contents.
            int32_t processDataDurMs = sdrParams.sysParams.coherentProcessingTimeMS   >
                                       sdrParams.sysParams.inCoherentProcessingTimeMS ?
                                       sdrParams.sysParams.coherentProcessingTimeMS   :
                                       sdrParams.sysParams.inCoherentProcessingTimeMS;


            double processDataDurSec = static_cast<double>(processDataDurMs) * 1e-3f;
            int32_t numSamplesPerCh = floor(processDataDurSec * samplingFreqHz);
            int32_t numSamplesPerFrame = numSamplesPerCh * channelSeperation;

            // Skip initial bytes
            fseek(fid, skipNumberOfBytes, SEEK_SET);

            // Read whole memory
            int32_t* p_RawData = new int32_t[numSamplesPerFrame];
            int32_t numWordsToRead = 0;
            int32_t numBitsPerSample = 2;
            int32_t numBitsPerSampleLog2 = log2(numBitsPerSample);
            int32_t numBitsPerSampleMask = (1 << numBitsPerSample) - 1;

            int32_t numBitsPerWord = (8 * sizeof(int32_t));
            int32_t numBitsPerWordLog2 = log2(numBitsPerWord);
            int32_t numSamplesPerWord = (numBitsPerWord >> numBitsPerSampleLog2);
            int32_t numSamplesPerWordLog2 = log2(numSamplesPerWord);
            int32_t numSamplesPerWordMask = (1 << numSamplesPerWordLog2) - 1;

            if (!dataType.compare("ubit2"))
            {
                numWordsToRead = (numSamplesPerFrame + 15) / numSamplesPerWord;
            }

            if (numWordsToRead)
            {
                size_t readFlag = fread(p_RawData, sizeof(int32_t), (unsigned int)numWordsToRead, fid);
            }
            else
            {
                ASSERT_WMSG("Number of samples to read from file must be non-zero");
            }

            if (!fileName.compare("Woodbine_47a"))
            {
                // TODO
            }
            else if (!fileName.compare("GPS_and_GIOVE_A-NN-fs16_3676-if4_1304.bin"))
            {
                // TODO
            }
            else if (!fileName.compare("GPSdata-DiscreteComponents-fs38_192-if9_55.bin"))
            {
                // TODO
            }
            else if (!fileName.compare("NTLab_Bands_GPS_GLONASS_L12.bin"))
            {
                for (int32_t chIdx = 0; chIdx < numRxDataChannels; chIdx++)
                {
                    // If memory is already allocated, re-use it.
                    if ((p_rxDataMem [chIdx].rxDataPerFrameMemType == PRE_PROCESS_MEM_ALLOCATED) &&
                        ((p_rxDataMem[chIdx].rxDataPerFrameBytesAllocated / sizeof(int32_t)) >= numSamplesPerCh)
                        )
                    {
                        // memory is sufficient, so use it.
                        // Re-use this memory but first reset it to zero.
                        memset(p_rxDataMem[chIdx].rxDataPerFrame, 0, sizeof(int32_t)* numSamplesPerCh);
                    }
                    else if ((p_rxDataMem[chIdx].rxDataPerFrameMemType == PRE_PROCESS_MEM_ALLOCATED) &&
                        ((p_rxDataMem[chIdx].rxDataPerFrameBytesAllocated / sizeof(int32_t)) < numSamplesPerCh)
                        )
                    {
                        delete [] p_rxDataMem[chIdx].rxDataPerFrame;

                        // Allocate memory if not already allocated.
                        p_rxDataMem[chIdx].rxDataPerFrame = new int32_t[numSamplesPerCh];
                        p_rxDataMem[chIdx].rxDataPerFrameMemType = PRE_PROCESS_MEM_ALLOCATED;
                        p_rxDataMem[chIdx].rxDataPerFrameBytesAllocated = (numSamplesPerCh * sizeof(int32_t));
                    }
                    else
                    {
                        // Allocate memory if not already allocated.
                        p_rxDataMem[chIdx].rxDataPerFrame = new int32_t[numSamplesPerCh];
                        p_rxDataMem[chIdx].rxDataPerFrameMemType = PRE_PROCESS_MEM_ALLOCATED;
                        p_rxDataMem[chIdx].rxDataPerFrameBytesAllocated = (numSamplesPerCh * sizeof(int32_t));
                    }
                    int32_t* p_data = p_rxDataMem[chIdx].rxDataPerFrame;
                    int32_t data;
                    int32_t sign;
                    int32_t val;
                    int32_t wordIdx;
                    int32_t bitShift;
                    int32_t bitIdx;
                    for (int32_t sampleIdx = 0; sampleIdx < numSamplesPerCh; sampleIdx++)
                    {
                        bitIdx = (chIdx + selectedChannel + sampleIdx * channelSeperation) << numBitsPerSampleLog2;
                        wordIdx = bitIdx >> numBitsPerWordLog2;
                        bitShift = (bitIdx & (numBitsPerWord-1));
                        data = (p_RawData[wordIdx] >> bitShift) & numBitsPerSampleMask;
                        
                        sign = 1 - ((data & 1) << 1);
                        val  = (((data >> 1) & 0x1) << 1)+ 1;

                        *p_data++ = sign * val;
                    }
                }
            }
            else
            {
                ASSERT_WMSG("Listed input data file is no supported by SDR. Add support in code!");
            }
            
            // Delete allocated data here
            delete[] p_RawData;
        }
    }
}
