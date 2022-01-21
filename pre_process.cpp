/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include "pre_process.h"
#include "pre_process_norm_acq_parcode.h"
#include "pre_process_weak_acq_dbzp.h"

#include "read_file_data.h"
#include "cfg/logger.h"
#include "cfg/ca_code_table.h"

namespace processing
{

	void PreProcess(
		SdrParams_t&         sdrParams,
		PreProcessSignals_t* p_prepSignal,
		int32_t              numAcqAlgos,
		RxDataChannelMem_t*  rxDataPerFrame,
		int32_t              numRxDataChannels
	)
	{

		// Read input data 
		ReadDataFile(sdrParams, rxDataPerFrame, numRxDataChannels);

		// Generate CA Code
		printString("s", "Generating C / A code.");

		// Instead use pre-defined one.
		const int8_t* caCodeTable = &CODE_TABLE[0][0];

		//  Iterate over algorithms to do
		// algorithm specific pre-processing

		for (int32_t algoIdx = 0; algoIdx < numAcqAlgos; algoIdx++)
		{
			if (!sdrParams.sysParams.acqAlgosList[algoIdx].compare("norm_acq_parcode"))
			{
				printString("s", "Preprocessing: pre_proc_norm_acq_parcode()");
				PreProcessNormAcqParcode(sdrParams, caCodeTable, &p_prepSignal[algoIdx]);
			}
			else if (!sdrParams.sysParams.acqAlgosList[algoIdx].compare("weak_acq_dbzp"))
			{
				printString("s", "Preprocessing: Calling pre_proc_weak_acq_dbzp()");
				PreProcessWeakAcqDbzp(sdrParams, caCodeTable, &p_prepSignal[algoIdx]);
			}
			else
			{
				ASSERT_WMSG("Acquisition algorithm is not recognized.");
			}


		}
		



	}


}
