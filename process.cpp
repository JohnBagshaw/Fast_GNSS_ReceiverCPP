/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include "process.h"
#include "process/acquisition/norm_acq_parcode.h"
#include "process/acquisition/weak_acq_dbzp.h"
#include "cfg/logger.h"

namespace processing
{

	void Process(
		SdrParams_t&         sdrParams,
		PreProcessSignals_t* p_prepSignal,
		int32_t              numAcqAlgos,
		RxDataChannelMem_t*  rxDataPerFrame,
		int32_t              numRxDataChannels,
		ProcessSignals_t*    procesSignals
	)
	{

		// Extract coherent frame data.
		
		int32_t currFileIdx     = sdrParams.stateParams.numFilesProcessed;
		int32_t numSamplesPerMs = int32_t(sdrParams.dataParamsList[currFileIdx].samplingFreqHz * 1e-3f);

		int32_t currFrameIdx        = sdrParams.stateParams.currFrameNum;
		int32_t numTotalFrames      = sdrParams.stateParams.numTotalFrames;
		int32_t  numSamplesPerFrame = sdrParams.sysParams.coherentProcessingTimeMS * numSamplesPerMs;
		int32_t  frameDataIdxStart  = currFrameIdx * numSamplesPerFrame;
		
		// currFrameNum * numSamplesPerFrame;
		string *acqAlgoList = sdrParams.sysParams.acqAlgosList;

		for (int32_t chIdx = 0; chIdx < numRxDataChannels; chIdx++)
		{

			RxDataChannelMem_t rxDataPerChannel = rxDataPerFrame[chIdx];
			rxDataPerChannel.rxDataPerFrame   = rxDataPerFrame[chIdx].rxDataPerFrame +
									              frameDataIdxStart;

			for (int32_t algoIdx = 0; algoIdx < numAcqAlgos; algoIdx++)
			{

				/** Print information to console
				*/

				string algoName = acqAlgoList[algoIdx];
				string infoStr; infoStr.reserve(200);
				sprintf(infoStr.data(),
					"Processing for frame :%d / %d, channel : %d / %d, algorithm : %s",
					sdrParams.stateParams.currFrameNum + 1,
					sdrParams.stateParams.numTotalFrames,
					chIdx+1,
					numRxDataChannels,
					algoName.c_str()
				);
				printString("s", infoStr.c_str());


				/** Invoke configured algorithms
				*/

				int32_t processResultIdx = chIdx * sdrParams.sysParams.numAcqAlgos + algoIdx;

				if (!algoName.compare("norm_acq_parcode"))
				{
					NormAcqParcode(
						sdrParams, 
						&p_prepSignal[algoIdx], 
						&rxDataPerChannel, 
						&procesSignals[processResultIdx]
					);
				}
				else if (!algoName.compare("weak_acq_dbzp"))
				{
					WeakAcqDbzp(
						sdrParams,
						&p_prepSignal[algoIdx],
						&rxDataPerChannel,
						&procesSignals[processResultIdx]
					);
				}
				else
				{
					ASSERT_WMSG("Acquisition algorithm is not recognized.");
				}
			}
		}
	}
}
