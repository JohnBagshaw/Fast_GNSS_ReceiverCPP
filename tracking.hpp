/**
 * @file tracking.hpp
 *
 * @brief C++ header file for SDR tracking
 *
 * Project Title: GNSS-R SDR
 * Author       : John Bagshaw
 * Contact      : jotshaw@yorku.ca
 * Supervisor   : Prof.Sunil Bisnath
 * Institution  : York University, Canada.
 **/

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

void tracking(
  FILE* fid,
  const Channel_t* InChannels,
  const Settings_t* settings,
  TrackResults_t* trackResults,
  Channel_t* outChannels
);





