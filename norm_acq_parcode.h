 /**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#ifndef ___NORM_ACQ_PARCODE__
#define ___NORM_ACQ_PARCODE__

#include "norm_acq_parcode.h"
#include "cfg/config_sdr_params.h"

using namespace std;
using namespace config;

namespace processing
{

	void NormAcqParcode(
		SdrParams_t& sdrParams,
		PreProcessSignals_t* p_prepSignal,
		RxDataChannelMem_t*  rxDataPerFrame,
		ProcessSignals_t*    procesSignals
	);
}



#endif // __PRE_PROCESS_NORM_ACQ_PARCODE__
