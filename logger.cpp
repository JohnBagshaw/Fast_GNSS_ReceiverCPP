/**
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

#include <iostream>
#include <iomanip>
#include <assert.h>

#include "logger.h"


using namespace std;

namespace logger
{
	void printMsg(
		string logType,
		string msgString,
		const char* callerFunc,
		const char* callerFile,
		string swVersionInfo
	)
	{
		cout << logType << ": " << swVersionInfo << setw(20)
			<< " : " << callerFile << "." << callerFunc << "() :> <strong>%s</strong> \n";
		if (!logType.compare("ASSERT"))
		{
			assert(0);
		}
	}


    void printString(const char* str1...) // C-style "const char* fmt, ..." is also valid
    {
        va_list args;
        va_start(args, str1);

        while (*str1 != '\0') {
            if (*str1 == 's') {
                const char* str = va_arg(args, const char*);
                std::cout << str << " ";
            }
            else if (*str1 == 't') {
                string str = va_arg(args, string);
                std::cout << str << " ";
            }
            else if (*str1 == 'd') {
                // note automatic conversion to integral type
                int i = va_arg(args, int);
                std::cout << i << " ";
            }
            else if (*str1 == 'f') {
                double d = va_arg(args, double);
                std::cout << d << " ";
            }
            ++str1;
        }
        std::cout << endl;

        va_end(args);
    }

}

