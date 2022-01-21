/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include <algorithm>
#include "post_process_norm_acq_parcode.h"
#include "fftw3.h"
#include "cfg/ca_code_table.h"

namespace processing
{

    static void sortDescend(
        double* dataIn,
        int32_t* dataIdx,
        int32_t numElems
    )
    {
        double tmp;
        int32_t tmpIdx;
        for (int32_t i = 0; i < numElems - 1; i++)
        {
            for (int32_t j = i + 1; j < numElems; j++)
            {
                if (dataIn[i] < dataIn[j])
                {
                    tmp = dataIn[i];
                    dataIn[i] = dataIn[j];
                    dataIn[j] = tmp;

                    tmpIdx     = dataIdx[i];
                    dataIdx[i] = dataIdx[j];
                    dataIdx[j] = tmpIdx;
                }
            }
        }
    }

    void PostProcessNormAcqParcode(
        SdrParams_t& sdrParams,
        ProcessSignals_t* p_processSignal,
        RxDataChannelMem_t* rxDataPerFrame,
        PostProcessResults_t* postProcessResults
    )
    {
        // Parameters
        int32_t currFileIdx = sdrParams.stateParams.numFilesProcessed;
        DataParams_t *dataFileParams = &sdrParams.dataParamsList[currFileIdx];


        int32_t averFactor = p_processSignal->averFactor;
        int32_t numCodeSamples = p_processSignal->numCodeSamples;
        int32_t numDopplerSamples = p_processSignal->numDopplerBins;
        int32_t dopplerResHz = p_processSignal->dopplerResHz;

        int32_t* satelliteList = sdrParams.sysParams.acqSatelliteList;
        int32_t maxNumAcqSatellites = sdrParams.sysParams.numAcqSatellites;
        int32_t numPrns = sdrParams.sysParams.numCfgSatellites;

        // Local buffers
        int32_t bPeakLocAcq[numTotalSatellites];
        int32_t bPeakValueAcq[numTotalSatellites];
        double  bPeakMetricAcq[numTotalSatellites];
        int32_t bacqPrnList[numTotalSatellites];


        int32_t bPeakLocNacq[numTotalSatellites];
        int32_t bPeakValueNacq[numTotalSatellites];
        double  bPeakMetricNacq[numTotalSatellites];
        int32_t bnacqPrnList[numTotalSatellites];

        double   chipRateHz = sdrParams.sysParams.caCodeChipRateHz;
        int32_t numChipsPerMs = floor(chipRateHz * 1e-3);

        int32_t numCodeSamplesDs = (numCodeSamples / averFactor);
        int32_t numDdmSamplesPerPrn = numCodeSamplesDs * numDopplerSamples;


        /** Keep DDM Accumulator buffer here.
         * 
         * 
         */

        int32_t numSamplesCorr = numDdmSamplesPerPrn * numPrns;
        if (sdrParams.stateParams.currFrameNum == 0)
        {
            // Allocate memory and assign ddm processing buffer here
            postProcessResults->numSamplesPerPrn = numDdmSamplesPerPrn;
            if ((postProcessResults->accumDdmMemType == POST_PROCESS_MEM_ALLOCATED) &&
               ((postProcessResults->accumDdmBytesAllocated / sizeof(complex<double>)) >= numSamplesCorr)
                )
            {
                // Re-use this memory but first reset it to zero.
                fill(postProcessResults->accumDdm,
                    postProcessResults->accumDdm + numSamplesCorr,
                    complex<double>(0.0, 0.0));
            }
            else if ((postProcessResults->accumDdmMemType == POST_PROCESS_MEM_ALLOCATED) &&
                    ((postProcessResults->accumDdmBytesAllocated / sizeof(complex<double>)) < numSamplesCorr)
                )
            {
                delete[] postProcessResults->accumDdm;
                postProcessResults->accumDdm = new complex<double>[numSamplesCorr];
                postProcessResults->accumDdmMemType = POST_PROCESS_MEM_ALLOCATED;
                postProcessResults->accumDdmBytesAllocated = numSamplesCorr * sizeof(complex<double>);
            }
            else
            {
                // Allocate memory if not already allocated.
                postProcessResults->accumDdm = new complex<double>[numSamplesCorr];
                postProcessResults->accumDdmMemType = POST_PROCESS_MEM_ALLOCATED;
                postProcessResults->accumDdmBytesAllocated = numSamplesCorr * sizeof(complex<double>);
            }

            // Just copy the memory first time.
            memcpy(postProcessResults->accumDdm, p_processSignal->ddMap, numSamplesCorr * sizeof(complex<double>));
        }
        else
        {
            for (int32_t sIdx = 0; sIdx < numSamplesCorr; sIdx++)
            {
                postProcessResults->accumDdm[sIdx] += p_processSignal->ddMap[sIdx];
            }
        }

        // Check if a satellite should be acquired
        // if its peak metric is greater than threshold.
        int32_t acqCount = 0;
        int32_t nacqCount = 0;
        for (int32_t prnIdx = 0; prnIdx < numPrns; prnIdx++)
        {
            // Get relevant DDM Map
            complex<double>* p_ddmMapPrn = &postProcessResults->accumDdm[prnIdx * numDdmSamplesPerPrn];

            double peakValue = (0.0, 0.0);
            int32_t peakLoc = 0;
            for (int32_t sIdx = 0; sIdx < numDdmSamplesPerPrn; sIdx++)
            {
                if (p_ddmMapPrn[sIdx].real() > peakValue)
                {
                    peakValue = p_ddmMapPrn[sIdx].real();
                    peakLoc = sIdx;
                }
            }

            // Calculate peak metric for all satellite signals
            int32_t dopplerShiftInt = floor(peakLoc / numCodeSamplesDs);
            int32_t codeDelay = (peakLoc % numCodeSamplesDs);

            // Second max peak
            double secondPeakValue = (0.0, 0.0);
            for (int32_t sIdx = 0; sIdx < (numCodeSamplesDs - 19); sIdx++)
            {
                int32_t idx = ((sIdx + 10 + codeDelay ) % numCodeSamplesDs);
                if (p_ddmMapPrn[dopplerShiftInt * numCodeSamplesDs + idx].real() > secondPeakValue)
                {
                    secondPeakValue = p_ddmMapPrn[dopplerShiftInt * numCodeSamplesDs + idx].real();
                }
            }
            double peakMetric = peakValue / secondPeakValue;

            if (peakMetric > sdrParams.sysParams.acqThreshold)
            {
                bPeakLocAcq[acqCount] = peakLoc;
                bPeakValueAcq[acqCount] = peakValue;
                bPeakMetricAcq[acqCount] = peakMetric;
                bacqPrnList[acqCount] = satelliteList[prnIdx];
                acqCount++;
            }
            else
            {
                bPeakLocNacq[nacqCount] = peakLoc;
                bPeakValueNacq[nacqCount] = peakValue;
                bPeakMetricNacq[nacqCount] = peakMetric;
                bnacqPrnList[nacqCount] = satelliteList[prnIdx];
                nacqCount++;
            }
        }


        // Refine results for acquired PRNs
        int32_t acqPrnCount = acqCount;
        if (acqPrnCount)
        {

            // if acquired satellites are greater than required, choose strongest.

            if (acqPrnCount > maxNumAcqSatellites)
            {
                double  tbPeakMetricAcq [numTotalSatellites];
                int32_t dataIdx         [numTotalSatellites];
                for (int32_t i = 0; i < acqPrnCount; i++)
                {
                    tbPeakMetricAcq[i] = bPeakMetricAcq[i];
                    dataIdx[i] = bacqPrnList[i];
                }
                sortDescend(tbPeakMetricAcq, dataIdx, acqPrnCount);

                // Only retain maxNumAcqSatellites
                // Then process acquired prn
                for (int32_t i = maxNumAcqSatellites; i < acqPrnCount; i++)
                {
                    for (int32_t j = 0; j < acqCount; j++)
                    {
                        if (dataIdx[i] == bacqPrnList[j])
                        {
                            bnacqPrnList[nacqCount] = bacqPrnList[j];
                            bPeakMetricNacq[nacqCount] = bPeakMetricAcq[j];
                            break;
                        }
                    }
                }

                // First process acquired prn
                int32_t acqCountTmp = 0;
                for (int32_t i = 0; i < acqPrnCount; i++)
                {
                    bool skip = false;
                    for (int32_t j = 0; j < acqPrnCount-maxNumAcqSatellites; j++)
                    {
                        if (dataIdx[maxNumAcqSatellites + j] == bacqPrnList[i])
                        {
                            skip = true;
                            break;
                        }
                    }

                    if (!skip)
                    {
                        bacqPrnList    [acqCountTmp]  = bacqPrnList[i];
                        bPeakMetricAcq [acqCountTmp] = bPeakMetricAcq[i];
                        bPeakLocAcq    [acqCountTmp] = bPeakLocAcq[i];
                        acqCountTmp++;
                    }
                }
                acqCount = acqCountTmp;
            }

            // Post processing for each acquired satellite
            for (int32_t aprxIdx = 0; aprxIdx < acqCount; aprxIdx++)
            {
                // Get prn number and location
                int32_t prnNum = bacqPrnList[aprxIdx];
                int32_t peakLoc = bPeakLocAcq[aprxIdx];
                int32_t prnIdx = prnNum - 1;

                // Get doppler and code shift
                int32_t dopplerShiftInt = floor(peakLoc / numCodeSamplesDs);
                int32_t codeDelay       = (peakLoc % numCodeSamplesDs);

                complex<double>* p_ddmMapPrn = &postProcessResults->accumDdm[prnIdx * numDdmSamplesPerPrn];

                // Assigning buffers
                        double*  p_dopplerSpectrum = new double[numDopplerSamples];
                complex<double>* p_rxData          = new complex<double>[numCodeSamples];
                complex<double>* p_rxDataDemod     = new complex<double>[numCodeSamples];
                complex<double>* p_caCode          = new complex<double>[numCodeSamples];
                complex<double>* p_cplxExp         = new complex<double>[numCodeSamples];



                /**
                * Calculate Frequency Estimate
                */ 

                for (int32_t dbinIdx = 0; dbinIdx < numDopplerSamples; dbinIdx++)
                {
                    p_dopplerSpectrum[dbinIdx] = p_ddmMapPrn[dbinIdx * numCodeSamplesDs + codeDelay].real();
                }

                // Coarse estimate
                int32_t mOneIdx = (dopplerShiftInt-1) % numDopplerSamples;
                int32_t midIdx  = (dopplerShiftInt  ) % numDopplerSamples;
                int32_t pOneIdx = (dopplerShiftInt+1) % numDopplerSamples;

                double mOneVal = sqrt(p_dopplerSpectrum[mOneIdx]);
                double midVal  = sqrt(p_dopplerSpectrum[midIdx]);
                double pOneVal = sqrt(p_dopplerSpectrum[pOneIdx]);

                double doppleShiftError = 0.5 * (mOneVal - pOneVal) /
                                          (mOneVal - 2 * midVal + pOneVal);
                dopplerShiftInt = dopplerShiftInt - floor(numDopplerSamples / 2);
                double dopplerShiftHz = ((double)dopplerShiftInt + doppleShiftError) * dopplerResHz;
                double ifFreqEst = dataFileParams->intermFreqHz + dopplerShiftHz;


                dopplerShiftHz = 0;
                int32_t codeDelayError = 0;
                int32_t rxDataAdjIdx = codeDelay * averFactor;
                if (rxDataAdjIdx < (rxDataPerFrame->rxDataPerFrameBytesAllocated / sizeof(int32_t)))
                {
                    /**
                    * Refine Code Delay Estimate
                    *
                    *
                    */

                    // Multiply rx data with doppler sinusoidal signal
                    const double piMult2 = 2.0 * acos(-1);
                    const std::complex<double> i(0, 1);

                    for (int32_t sIdx = 0; sIdx < numCodeSamples; sIdx++)
                    {
                        complex<double> cplxExp = exp(i * piMult2 *
                            (ifFreqEst / dataFileParams->samplingFreqHz) *
                            (double)sIdx);

                        // demod rx data
                        p_rxDataDemod[sIdx] = 
                            static_cast<double>(rxDataPerFrame->rxDataPerFrame[rxDataAdjIdx + sIdx]) * cplxExp;

                        // Save data
                        p_rxData[sIdx].real(rxDataPerFrame->rxDataPerFrame[rxDataAdjIdx + sIdx]);
                        p_cplxExp[sIdx] = cplxExp;
                    }

                    //FILE *fid = fopen("norm_cplx.txt", "w");
                    //for (int32_t j = 0; j < numCodeSamples; j++)
                    //{
                    //        fprintf(fid, "%.30f, %.30f\n", 
                    //                p_rxDataDemod[j].real(), 
                    //                p_rxDataDemod[j].imag());
                    //}
                    //fclose(fid);

                    // Generate CA Code
                    // Upsamples CA Code to match length.
                    const int8_t* caCodeTable = &CODE_TABLE[prnIdx][0];

                    double fracIndex = 0;
                    int32_t numSamplesPerRowCaCode = 1023;
                    double upFactor = chipRateHz / dataFileParams->samplingFreqHz;
                    for (int32_t sampleIdx = 0; sampleIdx < numCodeSamples; sampleIdx++)
                    {
                        int32_t caCodeIdx = static_cast<int32_t>(fracIndex * upFactor);
                        p_caCode[sampleIdx].real(static_cast<double>(caCodeTable[caCodeIdx]));
                        fracIndex += 1;
                    }

                    // Do correlation
                    int32_t maxLag = averFactor;
                    complex<double>* p_corrOut = new complex<double>[2 * maxLag + 1];
                    // negative lags
                    int32_t lag = maxLag;
                    for (size_t i = 0; i < maxLag; i++)
                    {
                        lag = maxLag - i;
                        for (int32_t j = 0; j < numCodeSamples; j++)
                        {
                            p_corrOut[i] += (p_rxDataDemod[j] * conj(p_caCode[lag]));
                            lag++;
                            if (lag == numCodeSamples)
                            {
                                break;
                            }
                        }
                    }
                    // zero lag
                    for (int32_t j = 0; j < numCodeSamples; j++)
                    {
                        p_corrOut[maxLag] += (p_rxDataDemod[j] * conj(p_caCode[j]));
                    }
                    // positive lags
                    lag = maxLag;
                    for (size_t i = maxLag + 1; i < 2 * maxLag; i++)
                    {
                        lag = i - maxLag;
                        for (int32_t j = 0; j < numCodeSamples; j++)
                        {
                            p_corrOut[i] += (p_rxDataDemod[lag] * conj(p_caCode[j]));
                            lag++;
                            if (lag == numCodeSamples)
                            {
                                break;
                            }
                        }
                    }

                    // Find max
                    codeDelayError = 0;
                    double maxVal = 0.0;
                    for (int32_t i = 0; i < 2 * maxLag + 1; i++)
                    {
                        if (abs(p_corrOut[i]) > maxVal)
                        {
                            maxVal = abs(p_corrOut[i]);
                            codeDelayError = i - maxLag;
                        }
                    }


                    /**
                    * Refine Doppler Estimate
                    *
                    *
                    *
                    */
                    for (int32_t sIdx = 0; sIdx < numCodeSamples; sIdx++)
                    {
                        p_rxData[sIdx] = p_rxData[sIdx] * p_caCode[sIdx] * p_cplxExp[sIdx];
                    }

                    double acqDopplerBW = sdrParams.sysParams.acqDopplerBwKhz * 1e3;
                    int32_t fftNumPts = 2 * (1 << int32_t(ceil(log2(numCodeSamples))));
                    double deltaF = dataFileParams->samplingFreqHz / fftNumPts;
                    int32_t pbins = int32_t(ceil(0.5 * acqDopplerBW / deltaF));

                    int32_t fAxisStart = 0.5 * fftNumPts - pbins + 1;
                    int32_t fAxisLength = 2 * pbins + 1;

                    double* fftFreqBins = new double[fAxisLength];
                    double* fftOutMag   = new double[fAxisLength];

                    for (int32_t i = 0; i < fAxisLength; i++)
                    {
                        fftFreqBins[i] = double(fAxisStart + i - floor(0.5 * fftNumPts) - 1) * deltaF;
                    }
                    // Assign memories for correlation through FFT
                    complex<double>* placeHolderMem = new complex<double>[fftNumPts];

                    // Define FFT Plan
                    fftw_plan caCodeFftPlanFft = fftw_plan_dft_1d(
                        (int)fftNumPts,
                        reinterpret_cast<fftw_complex*>(placeHolderMem),
                        reinterpret_cast<fftw_complex*>(placeHolderMem),
                        FFTW_FORWARD,
                        FFTW_MEASURE
                    );
                    memcpy(placeHolderMem, p_rxData, numCodeSamples * sizeof(complex<double>));
                    fftw_execute(caCodeFftPlanFft);
                    fftw_destroy_plan(caCodeFftPlanFft);

                    for (int32_t i = 0; i < fAxisLength; i++)
                    {
                        int32_t idx = (i + fftNumPts - fAxisLength / 2) % fftNumPts;
                        fftOutMag[i] = abs(placeHolderMem[idx]);
                    }
                    maxVal = 0.0;
                    int32_t fftMaxIndex = 0;
                    for (int32_t i = 0; i < fAxisLength; i++)
                    {
                        if (maxVal < fftOutMag[i])
                        {
                            maxVal = fftOutMag[i];
                            fftMaxIndex = i;
                        }
                    }

                    mOneIdx = (fftMaxIndex - 1) % fAxisLength;
                    midIdx = (fftMaxIndex     ) % fAxisLength;
                    pOneIdx = (fftMaxIndex + 1) % fAxisLength;

                    double mOneVal = sqrt(fftOutMag[mOneIdx]);
                    double midVal = sqrt(fftOutMag[midIdx]);
                    double pOneVal = sqrt(fftOutMag[pOneIdx]);

                    double doppleShiftError1 = 0.5 * (mOneVal - pOneVal) /
                        (mOneVal - 2 * midVal + pOneVal);
                    dopplerShiftHz = fftFreqBins[fftMaxIndex] + doppleShiftError1 * deltaF;

                    delete[]p_corrOut;
                    delete[]fftFreqBins;
                    delete[]fftOutMag;
                    delete[]placeHolderMem;
                }

                codeDelay = codeDelayError + codeDelay * averFactor + 1;
                ifFreqEst = ifFreqEst - dopplerShiftHz;

                // Save results per PRN
                postProcessResults->prnResults[aprxIdx].dopplerShiftHz = dopplerShiftHz;
                postProcessResults->prnResults[aprxIdx].codeDelay      = codeDelay;
                postProcessResults->prnResults[aprxIdx].peakMetric     = bPeakMetricAcq[aprxIdx];
                postProcessResults->prnResults[aprxIdx].estIfFreqHz    = ifFreqEst;
                postProcessResults->prnResults[aprxIdx].satellitePrn   = prnNum;

                delete[] p_dopplerSpectrum;
                delete[] p_rxData         ;
                delete[] p_rxDataDemod    ;
                delete[] p_caCode         ;
                delete[] p_cplxExp        ;

            } // Acquired PRN processing

            postProcessResults->numAcqSatellites = acqCount;
        }
    }
}