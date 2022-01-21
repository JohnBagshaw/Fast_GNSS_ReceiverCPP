/**
* @file acquisition.hpp
*
* @brief C++ header file for acquisition functions
*
* Project Title : GNSS - R SDR
* Author :        John Bagshaw
* Contact :       jotshaw@yorku.ca
* Supervisors:    Prof.Sunil Bisnath
* Institution :   York University, Canada.
** /

#pragma once

 /***********************************
 * Includes
 ***********************************/
#include "structures.hpp"
#include <stdio.h>

 /***********************************
 * Defines
 ***********************************/

 /***********************************
 * Macros
 ***********************************/

 /***********************************
 * Static Variables
 ***********************************/

 /***********************************
 * Static Function Definitions
 ***********************************/

 /***********************************
 * Public Function Declarations
 ***********************************/

void demodulation(double* longSignal,const Settings_t* settings);


bool acquisition(const double* longSignal,const Settings_t* settings,AcqResults_t* acqResults);





