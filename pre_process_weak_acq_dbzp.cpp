/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include "pre_process_weak_acq_dbzp.h"
#include "fftw3.h"
#include <iostream>

namespace processing
{

	void PreProcessWeakAcqDbzp(
		SdrParams_t& sdrParams,
		const int8_t* caCodeTable,
		PreProcessSignals_t* p_prepSignal
	)
	{
		/* Prepare parameters */

		int32_t fileIdx       = sdrParams.stateParams.numFilesProcessed;
		double samplingFreqHz = sdrParams.dataParamsList[fileIdx].samplingFreqHz;
		double intermFreqHz   = sdrParams.dataParamsList[fileIdx].intermFreqHz;

		double chipRateHz      = sdrParams.sysParams.caCodeChipRateHz;
		double acqDopplerBwKhz = sdrParams.sysParams.acqDopplerBwKhz;
		int32_t cpiTimeMs      = sdrParams.sysParams.coherentProcessingTimeMS;

		int32_t samplesPerCpi  = static_cast<int32_t>(samplingFreqHz * cpiTimeMs * 1e-3f);
		int32_t samplesPerMs   = static_cast<int32_t>(samplingFreqHz * 1e-3f);
		int32_t numChipsPerMs  = static_cast<int32_t>(chipRateHz * 1e-3);
		int32_t numCaCodeFolds = cpiTimeMs;

		int32_t* p_prnList    = sdrParams.sysParams.acqSatelliteList;
		int32_t numSatellites = sdrParams.sysParams.numCfgSatellites;


		if (floor(samplesPerCpi) != samplesPerCpi)
		{
			ASSERT_WMSG("Number of samples per millisecond can''t have decimal point.");
		}


		/**
		* Prepare per algorithm per block common signals
		*/

		// Determine maximum extent to which averaging can be done
		int32_t averFactor    = 1;
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

		/**
		* Prepare per algorithm per block common signals
		*/

		// Find out integer number of samples per blockand doppler bins
		double numSamplesPerCpi = floor(samplesPerCpi) / floor(averFactor);

		// number of doppler bins should be atleast 2 / t_cpi
		int32_t minDopplerBins = round(acqDopplerBwKhz * cpiTimeMs);
		bool    minFound       = false;
		int32_t numDopplerBins = INT32_MAX;
		for (double searchDbinIdx = minDopplerBins; searchDbinIdx <= 2 * minDopplerBins; searchDbinIdx++)
		{
			if (floor(numSamplesPerCpi / searchDbinIdx) == (numSamplesPerCpi / searchDbinIdx))
			{
				if (searchDbinIdx <= numDopplerBins)
				{
					minFound = true;
					numDopplerBins = static_cast<int32_t>(searchDbinIdx);
				}
			}
		}
		if (!minFound)
		{
			ASSERT_WMSG("Unable to find integer multiple of blockSize for DBZP");
		}

		int32_t numSamplesPerRowCaCode = 1023;
		int32_t numSamplesPerBlock     = static_cast<int32_t>(numSamplesPerCpi) / numDopplerBins;
		int32_t samplesPerMsDs         = static_cast<int32_t>(samplesPerMs / averFactor);
		double   upFactor              = numChipsPerMs * averFactor / (double)samplesPerMs;
		int32_t numBlocks              = samplesPerCpi / averFactor / numSamplesPerBlock;

		/**
		* Allocate memory for CA Code table in DBZP format.
		*/
		int32_t numCaCodeSamples = numSatellites * numBlocks * (2 * numSamplesPerBlock);
		if ((p_prepSignal->caCodeTableMemType == PRE_PROCESS_MEM_ALLOCATED) &&
			((p_prepSignal->caCodeTableBytesAllocated / sizeof(complex<double>)) >= numCaCodeSamples)
			)
		{
			// Re-use this memory but first reset it to zero.
			fill(p_prepSignal->dopplerCplxExp, 
				 p_prepSignal->dopplerCplxExp + numCaCodeSamples, 
				 complex<double>(0.0, 0.0));
		}
		else if ((p_prepSignal->caCodeTableMemType == PRE_PROCESS_MEM_ALLOCATED) &&
			    ((p_prepSignal->caCodeTableBytesAllocated / sizeof(complex<double>)) < numCaCodeSamples)
			)
		{
			delete[] p_prepSignal->caCodeTable;
			p_prepSignal->caCodeTable = new complex<double>[numCaCodeSamples];
			p_prepSignal->caCodeTableMemType = PRE_PROCESS_MEM_ALLOCATED;
			p_prepSignal->caCodeTableBytesAllocated = numCaCodeSamples * sizeof(complex<double>);
		}
		else
		{
			// Allocate memory if not already allocated.
			p_prepSignal->caCodeTable               = new complex<double>[numCaCodeSamples];
			p_prepSignal->caCodeTableMemType        = PRE_PROCESS_MEM_ALLOCATED;
			p_prepSignal->caCodeTableBytesAllocated = numCaCodeSamples * sizeof(complex<double>);
		}

		const int8_t* inCaCode = caCodeTable;
		complex<double>* outCaCode = p_prepSignal->caCodeTable;

		// Re-arrange in DBZP format.
		int32_t prnNum = 0;
		double fracIndex = 0;
		for (int32_t prnIdx = 0; prnIdx < numSatellites; prnIdx++)
		{
			// Upsamples CA Code to match length.
			prnNum = p_prnList[prnIdx] - 1;
			inCaCode = &caCodeTable[prnNum * numSamplesPerRowCaCode];

			fracIndex = 0;
			for (int32_t blockIdx = 0; blockIdx < numBlocks; blockIdx++)
			{
				for (int32_t sampleIdx = 0; sampleIdx < numSamplesPerBlock; sampleIdx++)
				{ 
					int32_t caCodeIdx = static_cast<int32_t>(upFactor * fracIndex) % numSamplesPerRowCaCode;
					outCaCode->real(static_cast<double>(inCaCode[caCodeIdx]));
					outCaCode += 1;
					fracIndex += 1;
				}
				outCaCode += numSamplesPerBlock;
			}
		}


		/**
		* Prepare complex exponential
		*/

		const double piMult2 = 2.0 * acos(-1);
		const std::complex<double> i(0, 1);

		//int32_t numDopplerSamples = floor(acqDopplerBwKhz * 1e3 / acqDopplerResHz) + 1;
		int32_t numDopplerSamples = 0;
		intermFreqHz;

		if ((p_prepSignal->dopplerCplxExpMemType == PRE_PROCESS_MEM_ALLOCATED) &&
			((p_prepSignal->dopplerCplxExpBytesAllocated / sizeof(complex<double>)) >= samplesPerCpi)
			)
		{
			// Re-use this memory but first reset it to zero.
			fill(p_prepSignal->dopplerCplxExp,
				p_prepSignal->dopplerCplxExp + samplesPerCpi,
				complex<double>(0.0, 0.0));
		}
		else if ((p_prepSignal->dopplerCplxExpMemType == PRE_PROCESS_MEM_ALLOCATED) &&
			((p_prepSignal->dopplerCplxExpBytesAllocated / sizeof(complex<double>)) < samplesPerCpi)
			)
		{
			delete[] p_prepSignal->dopplerCplxExp;
			p_prepSignal->dopplerCplxExp = new complex<double>[samplesPerCpi];
			p_prepSignal->dopplerCplxExpMemType = PRE_PROCESS_MEM_ALLOCATED;
			p_prepSignal->dopplerCplxExpBytesAllocated = samplesPerCpi * sizeof(complex<double>);

		}
		else
		{
			// Allocate memory if not already allocated.
			p_prepSignal->dopplerCplxExp = new complex<double>[samplesPerCpi];
			p_prepSignal->dopplerCplxExpMemType = PRE_PROCESS_MEM_ALLOCATED;
			p_prepSignal->dopplerCplxExpBytesAllocated = samplesPerCpi * sizeof(complex<double>);
		}

		complex<double>* outBuff = p_prepSignal->dopplerCplxExp;
		for (int32_t sampleIdx = 0; sampleIdx < samplesPerCpi; sampleIdx++)
		{
			*outBuff++ = exp(i * piMult2 * (intermFreqHz / samplingFreqHz) * (double)sampleIdx);
		}

		/**
		* Assing values.
		*/
		p_prepSignal->numBlocks = numBlocks;
		p_prepSignal->numSamplesPerBlock = numSamplesPerBlock;
		p_prepSignal->averFactor = averFactor;
	}
}
