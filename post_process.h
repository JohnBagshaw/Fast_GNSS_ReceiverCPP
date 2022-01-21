/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#ifndef __POST_PROCESS__
#define __POST_PROCESS__

#include "cfg/config_sdr_params.h"

using namespace std;
using namespace config;

namespace processing
{

	void PostProcess(
		SdrParams_t&         sdrParams,
		ProcessSignals_t*    procesSignals,
		int32_t              numAcqAlgos,
		RxDataChannelMem_t*  rxDataPerFrame,
		int32_t              numRxDataChannels,
		PostProcessResults_t* PostProcessResults_s
	);


}



#endif // __POST_PROCESS__
