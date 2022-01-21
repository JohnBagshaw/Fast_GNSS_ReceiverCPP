/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#ifndef __LOGGER__
#define __LOGGER__

#include <string>
#include <cstdarg>

using namespace std;

namespace logger
{

#define PRINT_MSG(str) printMsg("PRINT", str, __func__, __FILE__, str)
#define ASSERT_WMSG(str) printMsg("ASSERT", str, __func__, __FILE__, str)

void printMsg(
	string logType,
	string msgString,
	const char* callerFunc,
	const char* callerFile,
	string swVersionInfo
);

void printString(const char* str1...);

}



#endif // __LOGGER_
