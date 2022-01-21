/**
 * @file gen_code.hpp
 *
 * @brief code generation header file for SDR functions in C++
 *
 * Project Title: GNSS-R SDR
 * Author :       John Bagshaw
 * Contact :      jotshaw@yorku.ca
 * Supervisors:   Prof.Sunil Bisnath
 * Institution :  York University, Canada.
 **/

#pragma once

 /***********************************
 * Includes
 ***********************************/
#include "structures.hpp"

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

void generateCAcode(const Settings_t* settings,double* caCode, int PRN);

void makeCaTable(const Settings_t* settings, double**caCodesTable,int rows, int cols,double tS,double tC);