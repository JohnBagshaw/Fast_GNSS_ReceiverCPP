/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#ifndef __PRE_PROCESS_NORM_ACQ_PARCODE__
#define __PRE_PROCESS_NORM_ACQ_PARCODE__

#include "cfg/config_sdr_params.h"

using namespace std;
using namespace config;

namespace processing
{

	void PreProcessNormAcqParcode(
		SdrParams_t& sdrParams,
		const int8_t* caCodeTable,
		PreProcessSignals_t* p_prepSignal
	);


}



#endif // __PRE_PROCESS_NORM_ACQ_PARCODE__
