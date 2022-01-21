/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include "norm_acq_parcode.h"
#include "fftw3.h"
#include <algorithm> 
#include <iostream>

namespace processing
{

    void NormAcqParcode(
        SdrParams_t&         sdrParams,
        PreProcessSignals_t* p_prepSignal,
        RxDataChannelMem_t*  rxDataPerFrame,
        ProcessSignals_t*    procesSignals
    )
    {
        int32_t optimizeOption = p_prepSignal->optimizeOption;


        if (optimizeOption == 1)
        {
            // Parameters
            int32_t averFactor    = p_prepSignal->averFactor;
            int32_t* p_prnList    = sdrParams.sysParams.acqSatelliteList;
            int32_t numSatellites = sdrParams.sysParams.numCfgSatellites;

            int32_t numBlocks         = sdrParams.sysParams.coherentProcessingTimeMS;
            int32_t currFileIdx       = sdrParams.stateParams.numFilesProcessed;
            int32_t numCodeSamples    = sdrParams.dataParamsList[currFileIdx].samplingFreqHz * 1e-3;
            int32_t numDopplerSamples = floor(sdrParams.sysParams.acqDopplerBwKhz * 1e3 /
                                              sdrParams.sysParams.acqDopplerResHz) + 1;

            int32_t numSamplesMsDs = numCodeSamples / averFactor;
            int32_t rxDataStartIdx = 0;

            int32_t numSamplesCorr = numSatellites * numDopplerSamples * numSamplesMsDs;
            if (( procesSignals->ddMapMemType == PROCESS_MEM_ALLOCATED) &&
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
            complex<double>* demodRxDataCorr = procesSignals->ddMap;
            complex<double>* demodRxData = new complex<double>[numDopplerSamples * numSamplesMsDs];

            // Iterate over each block to do the processing
            for (int32_t blIdx = 0; blIdx < numBlocks; blIdx++)
            {
                // Input data placeholder
                complex<double>* p_demodRxData = demodRxData;

                // Multiply carrier signal with block data.
                // and take average...
                double multInvAverFactor    = 1.0 / double(averFactor);
                int32_t* rxDataPerBlock     = &rxDataPerFrame->rxDataPerFrame[rxDataStartIdx];
                complex<double>* dCplxExp    = p_prepSignal->dopplerCplxExp;
                complex<double> accum        = (0.0, 0.0);
                int32_t averFacCount        = 0;

                for (int32_t sampleIdx      = 0; sampleIdx < numCodeSamples * numDopplerSamples; sampleIdx++)
                {
                    int32_t sampleModSizeIdx = sampleIdx % numCodeSamples;
                    accum += dCplxExp[sampleIdx] * static_cast<double>(rxDataPerBlock[sampleModSizeIdx]);
                    averFacCount++;

                    if (averFacCount == averFactor)
                    {
                        *p_demodRxData++ = accum * multInvAverFactor;
                        averFacCount = 0;
                        accum.real(0.0);
                        accum.imag(0.0);
                    }
                }

                // Assign memories for correlation through FFT
                complex<double>* placeHolderMem = new complex<double>[numSamplesMsDs];

                // Define FFT Plan
                fftw_plan caCodeFftPlanFft = fftw_plan_dft_1d(
                    (int)numSamplesMsDs,
                    reinterpret_cast<fftw_complex*>(placeHolderMem),
                    reinterpret_cast<fftw_complex*>(placeHolderMem),
                    FFTW_FORWARD,
                    FFTW_MEASURE
                );

                fftw_plan caCodeFftPlanIfft = fftw_plan_dft_1d(
                    (int)numSamplesMsDs,
                    reinterpret_cast<fftw_complex*>(placeHolderMem),
                    reinterpret_cast<fftw_complex*>(placeHolderMem),
                    FFTW_BACKWARD,
                    FFTW_MEASURE
                );



                int32_t prnNum = 0;
                complex<double>* p_caCade = p_prepSignal->caCodeTable;
                double fftLenScale = 1.0 / double(numSamplesMsDs);
                complex<double>* p_demodRxDataCorr = demodRxDataCorr;

                for (int32_t prnIdx = 0; prnIdx < numSatellites; prnIdx++)
                {
                    // Upsamples CA Code to match length.
                    prnNum = p_prnList[prnIdx] - 1;

                    p_demodRxData = demodRxData;
                    for (int32_t dBinIdx = 0; dBinIdx < numDopplerSamples; dBinIdx++)
                    {
                        // Copy data from RxData to placeholder memory
                        memcpy(placeHolderMem, p_demodRxData, sizeof(complex<double>) * numSamplesMsDs);
                        
                        // Take FFT
                        fftw_execute(caCodeFftPlanFft);
                        
                        // Multiply
                        for (int32_t sampleIdx = 0; sampleIdx < numSamplesMsDs; sampleIdx++)
                        {
                        	placeHolderMem[sampleIdx] *= conj(p_caCade[sampleIdx]) * fftLenScale;
                        }
                        
                        // Take IFFT
                        fftw_execute(caCodeFftPlanIfft);

                        if (!blIdx)
                        {
                            // Add to previous block result. 
                            for (int32_t sampleIdx = 0; sampleIdx < numSamplesMsDs; sampleIdx++)
                            {
                                p_demodRxDataCorr[sampleIdx] = placeHolderMem[sampleIdx] * conj(placeHolderMem[sampleIdx]);
                            }
                        }
                        else
                        {
                            // Add to previous block result. 
                            for (int32_t sampleIdx = 0; sampleIdx < numSamplesMsDs; sampleIdx++)
                            {
                                p_demodRxDataCorr[sampleIdx] += placeHolderMem[sampleIdx] * conj(placeHolderMem[sampleIdx]);
                            }
                        }

                        // Move to next doppler bin
                        p_demodRxData += numSamplesMsDs;
                        p_demodRxDataCorr += numSamplesMsDs;
                    }

                    // Move to CA code for next PRN
                    p_caCade += numSamplesMsDs;
                }
                // Destroy FFT plans
                fftw_destroy_plan(caCodeFftPlanFft);
                fftw_destroy_plan(caCodeFftPlanIfft);
                delete[] placeHolderMem;

                // Next set of input data.
                rxDataStartIdx += numCodeSamples;
            }
            delete[]demodRxData;

            procesSignals->ddMap        = demodRxDataCorr;
            procesSignals->averFactor   = averFactor;
            procesSignals->dopplerResHz = p_prepSignal->dopplerResHz;
            procesSignals->numCodeSamples = numCodeSamples;
            procesSignals->numDopplerBins = numDopplerSamples;

		}
        else if (optimizeOption == 2)
		{
            // Parameters
            int32_t averFactor = p_prepSignal->averFactor;
            int32_t* p_prnList = sdrParams.sysParams.acqSatelliteList;
            int32_t numSatellites = sdrParams.sysParams.numCfgSatellites;
            
            int32_t numBlocks = sdrParams.sysParams.coherentProcessingTimeMS;
            int32_t currFileIdx = sdrParams.stateParams.numFilesProcessed;
            int32_t numCodeSamples = sdrParams.dataParamsList[currFileIdx].samplingFreqHz * 1e-3;
                  int32_t numDopplerSamples = p_prepSignal->numDopplerSamples;
            
            int32_t numSamplesMsDs = numCodeSamples / averFactor;
            int32_t rxDataStartIdx = 0;

            int32_t numSamplesCorr = numSatellites * numDopplerSamples * numSamplesMsDs;
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
            complex<double>* demodRxDataCorr = procesSignals->ddMap;
            complex<double>* demodRxData = new complex<double>[numSamplesMsDs];
            
            // Iterate over each block to do the processing
            for (int32_t blIdx = 0; blIdx < numBlocks; blIdx++)
			{
				// Input data placeholder
				complex<double>* p_demodRxData = demodRxData;

				// Multiply carrier signal with block data.
				// and take average...
				double multInvAverFactor = 1.0 / double(averFactor);
				int32_t* rxDataPerBlock = &rxDataPerFrame->rxDataPerFrame[rxDataStartIdx];
				complex<double>* dCplxExpInit = p_prepSignal->dopplerCplxExp;
				complex<double> accum = (0.0, 0.0);
				int32_t averFacCount = 0;

				for (int32_t sampleIdx = 0; sampleIdx < numSamplesMsDs; sampleIdx++)
				{
					accum.real(0.0);
					accum.imag(0.0);
					for (int32_t avIdx = 0; avIdx < averFactor; avIdx++)
					{
						accum += (dCplxExpInit[avIdx] * static_cast<double>(rxDataPerBlock[avIdx]));
					}
					rxDataPerBlock += averFactor;
					dCplxExpInit += averFactor;
					*p_demodRxData++ = accum * multInvAverFactor;
				}

				// Initial FFT

				// Assign memories for correlation through FFT
				complex<double>* placeHolderMem = new complex<double>[numSamplesMsDs];

				// Define FFT Plan
				fftw_plan caCodeFftPlanFft = fftw_plan_dft_1d(
					(int)numSamplesMsDs,
					reinterpret_cast<fftw_complex*>(placeHolderMem),
					reinterpret_cast<fftw_complex*>(placeHolderMem),
					FFTW_FORWARD,
					FFTW_ESTIMATE
				);
				p_demodRxData = demodRxData;
				memcpy(placeHolderMem, p_demodRxData, sizeof(complex<double>) * numSamplesMsDs);
				fftw_execute(caCodeFftPlanFft); // Take FFT
				fftw_destroy_plan(caCodeFftPlanFft);
				memcpy(p_demodRxData, placeHolderMem, sizeof(complex<double>) * numSamplesMsDs);

				// PRN processing
				fftw_plan caCodeFftPlanIfft = fftw_plan_dft_1d(
					(int)numSamplesMsDs,
					reinterpret_cast<fftw_complex*>(placeHolderMem),
					reinterpret_cast<fftw_complex*>(placeHolderMem),
					FFTW_BACKWARD,
					FFTW_ESTIMATE
				);

				int32_t prnNum = 0;
				complex<double>* p_caCade = p_prepSignal->caCodeTable;
				double fftLenScale = double(1.0) / double(numSamplesMsDs);
				complex<double>* p_demodRxDataCorr = demodRxDataCorr;
				int32_t shiftFactor = p_prepSignal->shiftFactor;

				for (int32_t prnIdx = 0; prnIdx < numSatellites; prnIdx++)
				{
					// Upsamples CA Code to match length.
					prnNum = p_prnList[prnIdx] - 1;

					p_demodRxData = demodRxData;
					for (int32_t dBinIdx = 0; dBinIdx < numDopplerSamples; dBinIdx++)
					{
						// Multiply
						int32_t shift = floor(shiftFactor * dBinIdx);
						int32_t rotIdxInit = (numSamplesMsDs - shift) % numSamplesMsDs;
						for (int32_t sampleIdx = 0; sampleIdx < numSamplesMsDs; sampleIdx++)
						{
							int32_t rotIdx = (rotIdxInit + sampleIdx) % numSamplesMsDs;
							placeHolderMem[sampleIdx] = p_demodRxData[rotIdx] * conj(p_caCade[sampleIdx]);
							placeHolderMem[sampleIdx] *= fftLenScale;
						}

						// Take IFFT
						fftw_execute(caCodeFftPlanIfft);

						if (!blIdx)
						{
							// Add to previous block result. 
							for (int32_t sampleIdx = 0; sampleIdx < numSamplesMsDs; sampleIdx++)
							{
								p_demodRxDataCorr[sampleIdx] = placeHolderMem[sampleIdx] * conj(placeHolderMem[sampleIdx]);
							}
						}
						else
						{
							// Add to previous block result. 
							for (int32_t sampleIdx = 0; sampleIdx < numSamplesMsDs; sampleIdx++)
							{
								p_demodRxDataCorr[sampleIdx] += placeHolderMem[sampleIdx] * conj(placeHolderMem[sampleIdx]);
							}
						}
						p_demodRxDataCorr += numSamplesMsDs;
					}

					// Move to CA code for next PRN
					p_caCade += numSamplesMsDs;
				}
				// Destroy FFT plans
				fftw_destroy_plan(caCodeFftPlanIfft);
				delete[] placeHolderMem;

				// Next set of input data.
				rxDataStartIdx += numCodeSamples;
			}
            delete[] demodRxData;

            procesSignals->averFactor = averFactor;
            procesSignals->dopplerResHz = p_prepSignal->dopplerResHz;
            procesSignals->numCodeSamples = numCodeSamples;
            procesSignals->numDopplerBins = numDopplerSamples;

        }
        else
        {
            ASSERT_WMSG("Normal parallel code search has invalid optimization option.");
        }
    }
}
