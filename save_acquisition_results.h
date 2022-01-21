
/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#ifndef __SAVE_ACQUISITION_RESULTS__
#define __SAVE_ACQUISITION_RESULTS__

#include "cfg/config_sdr_params.h"

using namespace std;
using namespace config;

namespace processing
{

	void SaveAcqResults(
		SdrParams_t&          sdrParams,
		int32_t               numAcqAlgos,
		int32_t               numRxDataChannels,
		PostProcessResults_t* PostProcessResults_s
	);


}



#endif // __POST_PROCESS__
