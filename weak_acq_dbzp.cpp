/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include "weak_acq_dbzp.h"
#include "cfg/logger.h"
#include "fftw3.h"

namespace processing
{
    void WeakAcqDbzp(
        SdrParams_t& sdrParams,
        PreProcessSignals_t* p_prepSignal,
        RxDataChannelMem_t* rxDataPerFrame,
        ProcessSignals_t* procesSignals
    )
    {

        complex<double> *p_caCodesTable = p_prepSignal->caCodeTable;
        complex<double>* p_dopplerFreqExp = p_prepSignal->dopplerCplxExp;

        int32_t numBlocks          = p_prepSignal->numBlocks;
        int32_t averFactor         = p_prepSignal->averFactor;
        int32_t numSamplesPerBlock = p_prepSignal->numSamplesPerBlock;

        int32_t* p_prnList    = sdrParams.sysParams.acqSatelliteList;
        int32_t numSatellites = sdrParams.sysParams.numCfgSatellites;

        int32_t currFileIdx    = sdrParams.stateParams.numFilesProcessed;
        double   samplingFreqHz = sdrParams.dataParamsList[currFileIdx].samplingFreqHz;

        int32_t numDopplerSamples = numBlocks;
        int32_t numCohIntMs       = sdrParams.sysParams.coherentProcessingTimeMS;

        int32_t dopplerDftInterpFactor = 4;
        int32_t numDopplerFftBins = dopplerDftInterpFactor * numDopplerSamples;
        int32_t numSamplesDsPerCpi = numSamplesPerBlock * numDopplerSamples;

        // Multiply rx data with complex exponentital to bring to baseband
        // and take average...
        
        double multInvAverFactor      = 1.0 / double(averFactor);
        int32_t* rxDataPerBlock      = rxDataPerFrame->rxDataPerFrame;
        complex<double>* dCplxExpInit = p_prepSignal->dopplerCplxExp;

        complex<double>* p_rxData = new complex<double>[2 * numSamplesDsPerCpi];
        complex<double> accum = (0.0, 0.0);
        complex<double>* rxData = p_rxData;
        for (int32_t blkIdx = 0; blkIdx < numBlocks; blkIdx++)
        {
            for (int32_t sampleIdx = 0; sampleIdx < numSamplesPerBlock; sampleIdx++)
            {
                accum.real(0.0);
                accum.imag(0.0);
                for (int32_t avIdx = 0; avIdx < averFactor; avIdx++)
                {
                    accum += (dCplxExpInit[avIdx] * static_cast<double>(rxDataPerBlock[avIdx]));
                }
                rxDataPerBlock += averFactor;
                dCplxExpInit   += averFactor;
                *rxData++ = accum * multInvAverFactor;
            }

            if (blkIdx)
            {
                rxData -= numSamplesPerBlock;
                for (int32_t sampleIdx = 0; sampleIdx < numSamplesPerBlock; sampleIdx++)
                {
                    *(rxData - numSamplesPerBlock) = *rxData;
                    rxData++;
                }
            }
            rxData += numSamplesPerBlock;
        }
        complex<double>* rxDataStrt = p_rxData;
        for (int32_t sampleIdx = 0; sampleIdx < numSamplesPerBlock; sampleIdx++)
        {
            *(rxData - numSamplesPerBlock) = *rxDataStrt++;
            rxData++;
        }

        // Process each PRN
        // Do correlation
        // Memory is arranged as first satellite, then block permutation samples, then all doppler samples (blocks)
        // i.e. each prn is a matrix of total (numBlocks*numSamplesPerBlock) rows and (numSamplesPerBlock*dopplerDftInterpFactor) columns
        int32_t numSamplesCorr = numSatellites * numSamplesPerBlock * numBlocks * numBlocks * dopplerDftInterpFactor;
        
        if ((procesSignals->ddMapMemType == PROCESS_MEM_ALLOCATED) &&
            ((procesSignals->ddMapBytesAllocated / sizeof(complex<double>)) >= numSamplesCorr)
            )
        {
            // Re-use this memory but first reset it to zero.
            fill(procesSignals->ddMap,
                procesSignals->ddMap + numSamplesCorr,
                complex<double>(0.0, 0.0));
        }
        else if ((procesSignals->ddMapMemType == PROCESS_MEM_ALLOCATED) &&
            ((procesSignals->ddMapBytesAllocated / sizeof(complex<double>)) < numSamplesCorr)
            )
        {
            delete[] procesSignals->ddMap;
            procesSignals->ddMap = new complex<double>[numSamplesCorr];
            procesSignals->ddMapMemType = PROCESS_MEM_ALLOCATED;
            procesSignals->ddMapBytesAllocated = numSamplesCorr * sizeof(complex<double>);
        }
        else
        {
            // Allocate memory if not already allocated.
            procesSignals->ddMap = new complex<double>[numSamplesCorr];
            procesSignals->ddMapMemType = PROCESS_MEM_ALLOCATED;
            procesSignals->ddMapBytesAllocated = numSamplesCorr * sizeof(complex<double>);
        }
        complex<double>* corrDdmMap = procesSignals->ddMap;


        // Memory is arranged as first satellite, then block, then all block permutation samples
        // i.e. each prn is a matrix of total (numBlock) rows and (numBlocks*numSamplesPerBlock) columns
        complex<double>* corrOut = new complex<double>[numSatellites * numSamplesPerBlock * numBlocks * numBlocks];
        complex<double>* caCodeFftPlaceHolderMem    = new complex<double>[2 * numSamplesPerBlock];
        complex<double>* rxDataFftPlaceHolderMem    = new complex<double>[2 * numSamplesPerBlock];
        complex<double>* corrDataIfftPlaceHolderMem = new complex<double>[2 * numSamplesPerBlock];
        complex<double>* corrDataDftPlaceHolderMem  = new complex<double>[numBlocks * dopplerDftInterpFactor];

        // Define FFT Plan
        fftw_plan caCodeFftPlanFft = fftw_plan_dft_1d(
            (int)(2 * numSamplesPerBlock),
            reinterpret_cast<fftw_complex*>(caCodeFftPlaceHolderMem),
            reinterpret_cast<fftw_complex*>(caCodeFftPlaceHolderMem),
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );
        fftw_plan rxDataFftPlanFft = fftw_plan_dft_1d(
            (int)(2 * numSamplesPerBlock),
            reinterpret_cast<fftw_complex*>(rxDataFftPlaceHolderMem),
            reinterpret_cast<fftw_complex*>(rxDataFftPlaceHolderMem),
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );
        fftw_plan corrDataIfftPlanFft = fftw_plan_dft_1d(
            (int)(2 * numSamplesPerBlock),
            reinterpret_cast<fftw_complex*>(corrDataIfftPlaceHolderMem),
            reinterpret_cast<fftw_complex*>(corrDataIfftPlaceHolderMem),
            FFTW_BACKWARD,
            FFTW_ESTIMATE
        );
        fftw_plan corrDataDftPlan = fftw_plan_dft_1d(
            (int)(numBlocks * dopplerDftInterpFactor),
            reinterpret_cast<fftw_complex*>(corrDataDftPlaceHolderMem),
            reinterpret_cast<fftw_complex*>(corrDataDftPlaceHolderMem),
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );

        // Alternate data buffer
        complex<double>* parCorrOut = new complex<double>[2*numSamplesPerBlock];
        complex<double>* p_parCorrOut = parCorrOut;

        int32_t halfDopplerSamples = numBlocks * dopplerDftInterpFactor / 2;
        complex<double>* rxDataFd  = new complex<double>[numBlocks * (2 * numSamplesPerBlock)];
        complex<double>* caCodeOut = new complex<double>[numBlocks * (2 * numSamplesPerBlock)];
        complex<double>* dftMem    = new complex<double>[2*halfDopplerSamples];

        // Re-arrange in DBZP format.
        int32_t prnNum = 0;
        double fracIndex = 0;
        double multInvFftLen = 1.0 / double(2 * numSamplesPerBlock);

        for (int32_t prnIdx = 0; prnIdx < numSatellites; prnIdx++)
        {
            complex<double>* p_corrOut   = &corrOut[prnIdx * numSamplesPerBlock * numBlocks * numBlocks];
            complex<double>* p_caCodeOut = caCodeOut;
            complex<double>* p_rxDataIn  = p_rxData;
            complex<double>* p_rxDataOut = rxDataFd;

            // FFT and IFFT
            complex<double>* p_caCodeIn = &p_caCodesTable[prnIdx * numBlocks * (2 * numSamplesPerBlock)];
            for (int32_t blkIdx = 0; blkIdx < numBlocks; blkIdx++)
            {
                memcpy(rxDataFftPlaceHolderMem, p_rxDataIn,   sizeof(complex<double>) * (2 * numSamplesPerBlock));
                memcpy(caCodeFftPlaceHolderMem, p_caCodeIn, sizeof(complex<double>)    * (2 * numSamplesPerBlock));
                fftw_execute(caCodeFftPlanFft); // Take IFFT
                fftw_execute(rxDataFftPlanFft); // Take FFT
                memcpy(p_caCodeOut, caCodeFftPlaceHolderMem, sizeof(complex<double>) * (2 * numSamplesPerBlock));
                memcpy(p_rxDataOut,  rxDataFftPlaceHolderMem, sizeof(complex<double>) * (2 * numSamplesPerBlock));

                p_caCodeOut += (2 * numSamplesPerBlock);
                p_rxDataOut += (2 * numSamplesPerBlock);

                p_rxDataIn += (2 * numSamplesPerBlock);
                p_caCodeIn  += (2 * numSamplesPerBlock);
            }

            complex<double>* p_corrOutBlk;
            for (int32_t permIdx = 0; permIdx < numBlocks; permIdx++)
            {

                // Get relevant addresses
                p_caCodeOut = caCodeOut;
                p_rxDataOut = rxDataFd;
                p_corrOutBlk = &p_corrOut[permIdx * numSamplesPerBlock];

                // Start index for PRN based offset to CA code.
                int32_t caCodeStartIdx = (numBlocks - permIdx) % numBlocks;

                // Iterate over the blocks
                for (int32_t blkIdx = 0; blkIdx < numBlocks; blkIdx++)
                {
                    int32_t caCodeBlockIdx = (caCodeStartIdx + blkIdx) % numBlocks;

                    // Permute the signal and multiply.
                    for (int32_t sIdx = 0; sIdx < (2 * numSamplesPerBlock); sIdx++)
                    {
                        p_parCorrOut[sIdx] = p_rxDataOut[sIdx] *
                                             conj(p_caCodeOut[caCodeBlockIdx * (2 * numSamplesPerBlock) + sIdx]) *
                                             multInvFftLen;
                    }
                    memcpy(corrDataIfftPlaceHolderMem, p_parCorrOut, sizeof(complex<double>) * (2 * numSamplesPerBlock));
                    fftw_execute(corrDataIfftPlanFft);
                    memcpy(p_corrOutBlk, corrDataIfftPlaceHolderMem, sizeof(complex<double>) * numSamplesPerBlock);
                    p_rxDataOut  += (2 * numSamplesPerBlock);
                    p_corrOutBlk += (numBlocks * numSamplesPerBlock);
                }
            }

            // Define memory to hold Doppler FFT output.
            complex<double>* p_corrDdmOut = &corrDdmMap[prnIdx * (numSamplesPerBlock * numBlocks * numBlocks * dopplerDftInterpFactor)];
            complex<double>* p_corrOutDplr = p_corrOut;


            // copy data from output to FFT input buffer
            // Iterate over each code bin (columns)
            int32_t numCodeBins = numSamplesPerBlock * numBlocks;
            for (int32_t codeIdx = 0; codeIdx < numCodeBins; codeIdx++)
            {
                // iterate over each row
                for (int32_t dBinIdx = 0; dBinIdx < numBlocks; dBinIdx++)
                {
                    int32_t dIdx = dBinIdx * numCodeBins + codeIdx;
                    p_corrDdmOut[dBinIdx] = p_corrOutDplr[dIdx];
                }

                // move to next doppler line
                p_corrDdmOut += (numBlocks * dopplerDftInterpFactor);
            }

            // OK Till now.


            //// Take FFT
            p_corrDdmOut = &corrDdmMap[prnIdx * (numSamplesPerBlock * numBlocks * numBlocks * dopplerDftInterpFactor)];
            for (int32_t codeIdx = 0; codeIdx < numCodeBins; codeIdx++)
            {
                memcpy(corrDataDftPlaceHolderMem, p_corrDdmOut, sizeof(complex<double>) * dopplerDftInterpFactor * numBlocks);
                fftw_execute(corrDataDftPlan);
                memcpy(p_corrDdmOut, corrDataDftPlaceHolderMem, sizeof(complex<double>) * dopplerDftInterpFactor* numBlocks);

                p_corrDdmOut += (dopplerDftInterpFactor * numBlocks);
            }

            // Take abs square
            // Copy to new buffer
            p_corrDdmOut = &corrDdmMap[prnIdx * (numSamplesPerBlock * numBlocks * numBlocks * dopplerDftInterpFactor)];

            // Correct till now.

            // Split buffer in half and add those two bufers
            int32_t halfCodeSamples = (numSamplesPerBlock * numBlocks * numBlocks * dopplerDftInterpFactor / 2);
            for (int32_t sampleIdx = 0; sampleIdx < halfCodeSamples; sampleIdx++)
            {
                 complex<double> data1 = p_corrDdmOut[sampleIdx] * conj(p_corrDdmOut[sampleIdx]);
                 complex<double> data2 = p_corrDdmOut[sampleIdx+ halfCodeSamples] * conj(p_corrDdmOut[sampleIdx+ halfCodeSamples]);

                 p_corrDdmOut[sampleIdx] = data1 + data2;
            }
            for (int32_t sampleIdx = halfCodeSamples; sampleIdx < 2 * halfCodeSamples ; sampleIdx++)
            {
                p_corrDdmOut[sampleIdx].real(0.0);
                p_corrDdmOut[sampleIdx].imag(0.0);
            }

            // FFT shift
            for (int32_t codeIdx = 0; codeIdx < (numSamplesPerBlock * numBlocks / 2); codeIdx++)
            {
                memcpy(&dftMem[halfDopplerSamples], p_corrDdmOut, sizeof(complex<double>) * halfDopplerSamples);
                memcpy(&dftMem[0], &p_corrDdmOut[halfDopplerSamples], sizeof(complex<double>)* halfDopplerSamples);
                memcpy(p_corrDdmOut, dftMem, sizeof(complex<double>)*2*halfDopplerSamples);

                p_corrDdmOut += (numBlocks * dopplerDftInterpFactor);
            }
        }

        // Define FFT Plan
        fftw_destroy_plan(caCodeFftPlanFft);
        fftw_destroy_plan(rxDataFftPlanFft);
        fftw_destroy_plan(corrDataIfftPlanFft);
        fftw_destroy_plan(corrDataDftPlan);
        delete[] parCorrOut;
        delete[] p_rxData;
        delete[] rxDataFd;
        delete[] caCodeOut;
        delete[] dftMem;
        delete[] corrOut;
        delete[] caCodeFftPlaceHolderMem;
        delete[] rxDataFftPlaceHolderMem;
        delete[] corrDataIfftPlaceHolderMem;
        delete[] corrDataDftPlaceHolderMem ;


        // Assign parameters
        procesSignals->ddMap = corrDdmMap;
        procesSignals->averFactor = averFactor;
        procesSignals->numSamplesPerBlock = numSamplesPerBlock;
        procesSignals->numBlocks = numBlocks / numCohIntMs;
        procesSignals->numDopplerBins = numDopplerSamples;
        procesSignals->dopplerResHz = (samplingFreqHz / averFactor / numSamplesPerBlock) / numDopplerFftBins;
        procesSignals->numDopplerFftBins = numDopplerFftBins;


     }
}
