/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#ifndef __POST_PROCESS_WEAK_ACQ_DBZP__
#define __POST_PROCESS_WEAK_ACQ_DBZP__

#include "cfg/config_sdr_params.h"

using namespace std;
using namespace config;

namespace processing
{

	void PostProcessWeakAcqDbzp(
		SdrParams_t&         sdrParams,
		ProcessSignals_t*     p_processSignal,
		RxDataChannelMem_t*   rxDataPerFrame,
		PostProcessResults_t* postProcessResults
	);


}



#endif // __POST_PROCESS_NORM_ACQ_PARCODE__
