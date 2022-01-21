/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include "cfg/logger.h"

#include "post_process.h"
#include "post_process_norm_acq_parcode.h"
#include "post_process_weak_acq_dbzp.h"

namespace processing
{

	void PostProcess(
		SdrParams_t&         sdrParams,
		ProcessSignals_t*    procesSignals,
		int32_t              numAcqAlgos,
		RxDataChannelMem_t*  rxDataPerFrame,
		int32_t              numRxDataChannels,
		PostProcessResults_t* PostProcessResults_s
	)
	{

		// Extract coherent frame data.
		int32_t currFileIdx = sdrParams.stateParams.numFilesProcessed;
		int32_t numChannels = sdrParams.dataParamsList[currFileIdx].totalChannels;
		int32_t selectedChannel = sdrParams.dataParamsList[currFileIdx].selectedChannel;


		// currFrameNum * numSamplesPerFrame;
		string* acqAlgoList = sdrParams.sysParams.acqAlgosList;

		for (int32_t chIdx = 0; chIdx < numRxDataChannels; chIdx++)
		{

			RxDataChannelMem_t rxDataPerChannel = rxDataPerFrame[chIdx];
			for (int32_t algoIdx = 0; algoIdx < numAcqAlgos; algoIdx++)
			{

				/** Print information to console
				*/

				string algoName = acqAlgoList[algoIdx];
				string infoStr; infoStr.reserve(200);
				sprintf(infoStr.data(),
					"Post Processing for frame :%d / %d, channel : %d / %d, algorithm : %s",
					sdrParams.stateParams.currFrameNum + 1,
					sdrParams.stateParams.numTotalFrames,
					chIdx + 1,
					numRxDataChannels,
					algoName.c_str()
				);
				printString("s", infoStr.c_str());


				/** Invoke configured algorithms
				*/

				int32_t processResultIdx = chIdx * sdrParams.sysParams.numAcqAlgos + algoIdx;

				if (!algoName.compare("norm_acq_parcode"))
				{
					PostProcessNormAcqParcode(
						sdrParams,
						&procesSignals[processResultIdx],
						&rxDataPerChannel,
						&PostProcessResults_s[processResultIdx]
					);
				}
				else if (!algoName.compare("weak_acq_dbzp"))
				{
 					PostProcessWeakAcqDbzp(
						sdrParams,
						&procesSignals[processResultIdx],
						&rxDataPerChannel,
						&PostProcessResults_s[processResultIdx]
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
