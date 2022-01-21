/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include "pre_process_norm_acq_parcode.h"
#include "fftw3.h"
#include <iostream>
#include <chrono>

namespace processing
{

    void PreProcessNormAcqParcode(
        SdrParams_t& sdrParams,
        const int8_t* caCodeTable,
        PreProcessSignals_t* p_prepSignal
    )
    {
        p_prepSignal->optimizeOption = 2; // 1,2 are only options

        /* Prepare parameters */

        int32_t fileIdx   = sdrParams.stateParams.numFilesProcessed;
		double chipRateHz = sdrParams.sysParams.caCodeChipRateHz;
		double acqDopplerBwKhz = sdrParams.sysParams.acqDopplerBwKhz;
		double acqDopplerResHz = sdrParams.sysParams.acqDopplerResHz;

		string fileName = sdrParams.stateParams.fileNames[fileIdx];
		double samplingFreqHz = sdrParams.dataParamsList[fileIdx].samplingFreqHz;
		double intermFreqHz = sdrParams.dataParamsList[fileIdx].intermFreqHz;


		int32_t samplesPerMs = static_cast<int32_t>(samplingFreqHz * 1e-3f);
		int32_t numChipsPerMs = static_cast<int32_t>(chipRateHz * 1e-3f);
		int32_t* p_prnList = sdrParams.sysParams.acqSatelliteList;
		int32_t numSatellites = sdrParams.sysParams.numCfgSatellites;

		if (floor(samplesPerMs) != samplesPerMs)
		{
			ASSERT_WMSG("Number of samples per millisecond can''t have decimal point.");
		}


		/**
		* Prepare per algorithm per block common signals
		*/

		// Determine maximum extent to which averaging can be done
		int32_t averFactor = 1;
		int32_t temAverFactor = 1;
		while ((samplingFreqHz / temAverFactor) > sdrParams.sysParams.minSamplingFreqHz)
		{
			temAverFactor++;
			double integerAverFactor = double(samplesPerMs) / temAverFactor;
			if (floor(integerAverFactor) == integerAverFactor)
			{
				averFactor = temAverFactor;
			}
		}

		// Allocate memory for FFT of CA Code
		int32_t samplesPerMsDs = static_cast<int32_t>(samplesPerMs / averFactor);
		int32_t numSamplesPerCaCode = samplesPerMsDs * numSatellites;

		// If memory is already allocated, re-use it.
		if ((p_prepSignal->caCodeTableMemType == PRE_PROCESS_MEM_ALLOCATED) &&
			((p_prepSignal->caCodeTableBytesAllocated / sizeof(complex<double>)) >= numSamplesPerCaCode)
			)
		{
			// Re-use this memory but first reset it to zero.
			fill(p_prepSignal->caCodeTable, p_prepSignal->caCodeTable+ numSamplesPerCaCode, complex<double>(0.0,0.0));
		}
		else if ((p_prepSignal->caCodeTableMemType == PRE_PROCESS_MEM_ALLOCATED) &&
			((p_prepSignal->caCodeTableBytesAllocated / sizeof(complex<double>)) < numSamplesPerCaCode)
			)
		{
			delete[] p_prepSignal->caCodeTable;
			p_prepSignal->caCodeTable = new complex<double>[numSamplesPerCaCode];
			p_prepSignal->caCodeTableMemType = PRE_PROCESS_MEM_ALLOCATED;
			p_prepSignal->caCodeTableBytesAllocated = numSamplesPerCaCode * sizeof(complex<double>);
		}
		else
		{
			// Allocate memory if not already allocated.
			p_prepSignal->caCodeTable = new complex<double>[numSamplesPerCaCode];
			p_prepSignal->caCodeTableMemType = PRE_PROCESS_MEM_ALLOCATED;
			p_prepSignal->caCodeTableBytesAllocated = numSamplesPerCaCode * sizeof(complex<double>);
		}

		/**
		* Prepare per algorithm per block common signals
		* each PRN C / A code sequence frequency transform here.
		*/

		int32_t numSamplesPerRowCaCode = 1023;
		double upFactor = chipRateHz * averFactor / samplingFreqHz;

		// Define FFT placeholder memory.
		complex<double>* placeHolderMem = new complex<double>[samplesPerMsDs];

		// cDefine FFT Plan
		fftw_plan caCodeFftPlan = fftw_plan_dft_1d(
			(int)samplesPerMsDs,
			reinterpret_cast<fftw_complex*>(placeHolderMem),
			reinterpret_cast<fftw_complex*>(placeHolderMem),
			FFTW_FORWARD, FFTW_ESTIMATE);

		const int8_t* inCaCode = caCodeTable;
		complex<double>* outCaCode = p_prepSignal->caCodeTable;

		int32_t prnNum = 0;
		double fracIndex = 0;
		for (int32_t prnIdx = 0; prnIdx < numSatellites; prnIdx++)
		{
			// Upsamples CA Code to match length.
			prnNum = p_prnList[prnIdx] - 1;
			inCaCode = &caCodeTable[prnNum * numSamplesPerRowCaCode];

			fracIndex = 0;
			for (int32_t sampleIdx = 0; sampleIdx < samplesPerMsDs; sampleIdx++)
			{
				int32_t caCodeIdx = static_cast<int32_t>(fracIndex * upFactor);
				outCaCode[sampleIdx].real(static_cast<double>(inCaCode[caCodeIdx]));
				fracIndex += 1;
			}

			// FFT.
			memcpy(placeHolderMem, outCaCode, sizeof(complex<double>) * samplesPerMsDs);
			fftw_execute(caCodeFftPlan);
			memcpy(outCaCode, placeHolderMem, sizeof(complex<double>) * samplesPerMsDs);

			// Move to next PRN
			outCaCode += samplesPerMsDs;
		}
		fftw_destroy_plan(caCodeFftPlan);
		delete[] placeHolderMem;


		const double piMult2 = 2.0 * acos(-1);
		const std::complex<double> i(0, 1);

		if (p_prepSignal->optimizeOption == 1)
		{
			int32_t numDopplerSamples = floor(acqDopplerBwKhz * 1e3 / acqDopplerResHz) + 1;
			intermFreqHz;

			p_prepSignal->dopplerCplxExp = new complex<double>[numDopplerSamples * samplesPerMs];

			complex<double>* outBuff = p_prepSignal->dopplerCplxExp;
			double dopplerFreq = -((acqDopplerBwKhz * 1000) * 0.5) + intermFreqHz;

			for (int32_t dBinIdx = 0; dBinIdx < numDopplerSamples; dBinIdx++)
			{
				for (int32_t sampleIdx = 0; sampleIdx < samplesPerMs; sampleIdx++)
				{
					*outBuff++ = exp(i * piMult2 * (dopplerFreq / samplingFreqHz) * (double)sampleIdx);
				}

				// Next doppler frequency
				dopplerFreq += acqDopplerResHz;
			}
			p_prepSignal->dopplerResHz = acqDopplerResHz;
		}
		else if (p_prepSignal->optimizeOption == 2)
		{
			int32_t acqDopplerResHz = 1000;
			int32_t numDopplerSamples = floor(acqDopplerBwKhz * 1e3 / acqDopplerResHz) + 1;
			double iFfreqVec = intermFreqHz - acqDopplerBwKhz * 0.5e3;

			/**
			* Allocate Complex Exponential Memories.
			*
			*/

			// If memory is already allocated, re-use it.
			if ((p_prepSignal->dopplerCplxExpMemType == PRE_PROCESS_MEM_ALLOCATED) &&
				((p_prepSignal->dopplerCplxExpBytesAllocated / sizeof(complex<double>)) >= samplesPerMs)
				)
			{
				// Re-use this memory but first reset it to zero.
				fill(p_prepSignal->dopplerCplxExp, p_prepSignal->dopplerCplxExp + samplesPerMs, complex<double>(0.0, 0.0));
			}
			else if ((p_prepSignal->dopplerCplxExpMemType == PRE_PROCESS_MEM_ALLOCATED) &&
				((p_prepSignal->dopplerCplxExpBytesAllocated / sizeof(complex<double>)) < samplesPerMs)
				)
			{
				delete[] p_prepSignal->dopplerCplxExp;
				p_prepSignal->dopplerCplxExp = new complex<double>[samplesPerMs];
				p_prepSignal->dopplerCplxExpMemType = PRE_PROCESS_MEM_ALLOCATED;
				p_prepSignal->dopplerCplxExpBytesAllocated = samplesPerMs * sizeof(complex<double>);
			}
			else
			{
				// Allocate memory if not already allocated.
				p_prepSignal->dopplerCplxExp = new complex<double>[samplesPerMs];
				p_prepSignal->dopplerCplxExpMemType = PRE_PROCESS_MEM_ALLOCATED;
				p_prepSignal->dopplerCplxExpBytesAllocated = samplesPerMs * sizeof(complex<double>);
			}

			// If memory is already allocated, re-use it.
			if ((p_prepSignal->dopplerCplxExpDeltaMemType == PRE_PROCESS_MEM_ALLOCATED) &&
				((p_prepSignal->dopplerCplxExpDeltaBytesAllocated / sizeof(complex<double>)) >= samplesPerMs)
				)
			{
				// Re-use this memory but first reset it to zero.
				fill(p_prepSignal->dopplerCplxExpDelta, p_prepSignal->dopplerCplxExpDelta + numSamplesPerCaCode, complex<double>(0.0, 0.0));
			}
			else if ((p_prepSignal->dopplerCplxExpDeltaMemType == PRE_PROCESS_MEM_ALLOCATED) &&
				((p_prepSignal->dopplerCplxExpDeltaBytesAllocated / sizeof(complex<double>)) < samplesPerMs)
				)
			{
				delete[] p_prepSignal->dopplerCplxExpDelta;
				p_prepSignal->dopplerCplxExpDelta = new complex<double>[samplesPerMs];
				p_prepSignal->dopplerCplxExpDeltaMemType = PRE_PROCESS_MEM_ALLOCATED;
				p_prepSignal->dopplerCplxExpDeltaBytesAllocated = samplesPerMs * sizeof(complex<double>);
			}
			else
			{
				// Allocate memory if not already allocated.
				p_prepSignal->dopplerCplxExpDelta = new complex<double>[samplesPerMs];
				p_prepSignal->dopplerCplxExpDeltaMemType = PRE_PROCESS_MEM_ALLOCATED;
				p_prepSignal->dopplerCplxExpDeltaBytesAllocated = samplesPerMs * sizeof(complex<double>);
			}

			complex<double>* outBuffInit = p_prepSignal->dopplerCplxExp;
			complex<double>* outBuffDelta = p_prepSignal->dopplerCplxExpDelta;

			for (int32_t sampleIdx = 0; sampleIdx < samplesPerMs; sampleIdx++)
			{
				*outBuffInit++ = exp(i * piMult2 * (iFfreqVec / samplingFreqHz) * (double)sampleIdx);
				*outBuffDelta++ = exp(i * piMult2 * (acqDopplerResHz / samplingFreqHz) * (double)sampleIdx);
			}

			// Assign internally used variables.
			p_prepSignal->shiftFactor = acqDopplerResHz * samplesPerMs / samplingFreqHz;
			p_prepSignal->numDopplerSamples = numDopplerSamples;
			p_prepSignal->dopplerResHz = acqDopplerResHz;

		}
		else
		{
			ASSERT_WMSG("Invalid optimization option");
		}

		// Assign average factor used.
		p_prepSignal->averFactor = averFactor;
	}
}
