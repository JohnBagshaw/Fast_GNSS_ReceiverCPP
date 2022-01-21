 /**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#ifndef __PRE_PROCESS_WEAK_ACQ_DBZP__
#define __PRE_PROCESS_WEAK_ACQ_DBZP__

#include "cfg/config_sdr_params.h"

using namespace std;
using namespace config;

namespace processing
{

void PreProcessWeakAcqDbzp(
	SdrParams_t& sdrParams,
	const int8_t* caCodeTable,
	PreProcessSignals_t* p_prepSignal
);

}

#endif // __PRE_PROCESS_WEAK_ACQ_DBZP__
