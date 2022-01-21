/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>

#include "cfg/logger.h"

#include "save_acquisition_results.h"

namespace processing
{

	void SaveAcqResults(
		SdrParams_t& sdrParams,
		int32_t               numAcqAlgos,
		int32_t               numRxDataChannels,
		PostProcessResults_t* PostProcessResults_s
	)
	{

		// Extract coherent frame data.
		int32_t currFrameIdx = sdrParams.stateParams.currFrameNum;
		int32_t currFileIdx = sdrParams.stateParams.numFilesProcessed;
		int32_t numChannels = sdrParams.dataParamsList[currFileIdx].totalChannels;
		int32_t selectedChannel = sdrParams.dataParamsList[currFileIdx].selectedChannel;
		string fileName = sdrParams.stateParams.fileNames[currFileIdx];
		string* acqAlgoList = sdrParams.sysParams.acqAlgosList;

		for (int32_t chIdx = 0; chIdx < numRxDataChannels; chIdx++)
		{
			for (int32_t algoIdx = 0; algoIdx < numAcqAlgos; algoIdx++)
			{
				// Extract relevan results structure
				int32_t pResultIdx = algoIdx + chIdx * numAcqAlgos;

				/** Print information to console
				*/
				string algoName = acqAlgoList[algoIdx];
				string infoStr; infoStr.reserve(200);
				sprintf(infoStr.data(), "Displaying Results======(File: %s, Channel: %d, Algorithm: %s)======",
					fileName.c_str(), chIdx+1, acqAlgoList[algoIdx].c_str());
				printString("s", infoStr.c_str());

				printString("s", "=========================================================================");
				sprintf(infoStr.data(), "Displaying Results======( %4d Satellites Detected )======",
					PostProcessResults_s[pResultIdx].numAcqSatellites);
				printString("s", infoStr.c_str());
				printString("s", "=========================================================================");


				std::ostringstream oss;
				auto t  = std::time(nullptr);
				auto tm = *std::localtime(&t);
				oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
				string timeStr = oss.str();

				size_t lastindex = fileName.find_last_of(".");
				string rawname = fileName.substr(0, lastindex);


				string outFileName; outFileName.reserve(500);
				sprintf(outFileName.data(), "file_%s_ch_%d_frame_%d_%s_%s.txt",
					rawname.c_str(),chIdx+1, currFrameIdx,acqAlgoList[algoIdx].c_str(),timeStr.c_str()
				);

				string fileOutFullName = sdrParams.stateParams.dataOutPath + string(outFileName.c_str());
				FILE* fid = fopen(fileOutFullName.c_str(), "w");
				if (fid == NULL)
				{
					ASSERT_WMSG("Output data file could not be opened.");
				}
				else
				{
					fprintf(fid, "NumAcqSatellites=%d\n", PostProcessResults_s[pResultIdx].numAcqSatellites);
					for (int32_t acqPrnIdx = 0; 
						acqPrnIdx < PostProcessResults_s[pResultIdx].numAcqSatellites; 
						acqPrnIdx++)
					{
						fprintf(fid, "{\n");
						fprintf(fid, "PRN            = %d\n",    PostProcessResults_s[pResultIdx].prnResults[acqPrnIdx].satellitePrn);
						fprintf(fid, "peakMetric     = %.15f\n", PostProcessResults_s[pResultIdx].prnResults[acqPrnIdx].peakMetric);
						fprintf(fid, "IfFreqHz       = %.15f\n", PostProcessResults_s[pResultIdx].prnResults[acqPrnIdx].estIfFreqHz);
						fprintf(fid, "codeDelay      = %d\n",    PostProcessResults_s[pResultIdx].prnResults[acqPrnIdx].codeDelay);
						fprintf(fid, "dopplerShiftHz = %.15f\n", PostProcessResults_s[pResultIdx].prnResults[acqPrnIdx].dopplerShiftHz);
						fprintf(fid, "}\n");

						string logStr; logStr.reserve(500);
						sprintf(logStr.data(), "(PRN: %4d, Metric: %10.4f, Code: %7d, Frequency: %10.4f)",
							PostProcessResults_s[pResultIdx].prnResults[acqPrnIdx].satellitePrn,
							PostProcessResults_s[pResultIdx].prnResults[acqPrnIdx].peakMetric,
							PostProcessResults_s[pResultIdx].prnResults[acqPrnIdx].codeDelay,
							PostProcessResults_s[pResultIdx].prnResults[acqPrnIdx].estIfFreqHz);
						printString("s", logStr.c_str());

					}
				}
				fclose(fid);
				printString("s", "=========================================================================");

			} // algo
		} // channel
	}
}
