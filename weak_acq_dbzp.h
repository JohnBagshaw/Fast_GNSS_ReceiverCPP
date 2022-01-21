/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#ifndef ___WEAK_ACQ_DBZP__
#define ___WEAK_ACQ_DBZP__

#include "cfg/config_sdr_params.h"

using namespace std;
using namespace config;

namespace processing
{
	void WeakAcqDbzp(
		SdrParams_t&         sdrParams,
		PreProcessSignals_t* p_prepSignal,
		RxDataChannelMem_t*  rxDataPerFrame,
		ProcessSignals_t*    procesSignals
	);
}



#endif // ___WEAK_ACQ_DBZP__
