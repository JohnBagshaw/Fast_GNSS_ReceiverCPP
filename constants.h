/**
* @file constants.h
*
* @brief C++ header file for SDR constants
*
* Project Title: GNSS-R SDR
* Author :       John Bagshaw
* Contact :      jotshaw@yorku.ca
* Supervisors:   Prof.Sunil Bisnath
* Institution :  York University, Canada.
**/


#pragma once

#define _USE_MATH_DEFINES

#define NUM_CHANNELS             (8)
#define MAX_FILE_NAME_LEN_CHARS (200)
#define C                       (2.99792458E+8)

#define Sig_ROUND(x) (x >= 0 ? 1 : -1)
#define GPS_SIN(x) (x > 0 ? 1 : -1)
#define GPS_COS(x) (x > M_PI/2 || x < -M_PI/2 ? -1 : 1)
//#define GPS_SINC(x) (x > 0x10000000 ? 0 : 1) //‭268435456‬
//#define GPS_COSC(x) (x >= 0x11000000 || x <= 0x00111111 ? 0 : 1)
#define GPS_SINC(y) ((y >= 0x080000000) ? 0 : 1) //‭80000000‬
#define GPS_COSC(y) (y >= 0x40000000 && y <= 0xc0000000 ? 0 : 1) // // //1073741824 //2147483678 // 1)1,073,741,824 2)2,147,483,648 3)3,221,225,472 4)4,294,967,296‬
#define UNWRAP_ANGLE(x) { if (x > M_PI) x -= 2 * M_PI; else if (x < -M_PI) x += 2 * M_PI;}

#define ACQ_MS        1
// These gains calculated "roughly" using the formulas provided in Paul Groves, Artech book
// tracking loop default gains
#define KCO_DEFAULT			0.004 // Kco = 4*BW*CODE_DURATION
#define KCO2_DEFAULT		1.0

// for the first half second, very high BW/Gain to pull the freq in fast
#define KCA2_FLL2_DEFAULT	1.15 //1.15	// with normalization, BW 50
// for the next half second, drop the loopbandwidth to better determine the freq before jump to PLL
#define KCA2_FLL1_DEFAULT	0.1//0.10	// with normalization, BW 15

#define KCA2_PLL_DEFAULT	0.10	// with normalization, BW 15
#define KCA3_PLL_DEFAULT	0.93	// with normalization, BW 15

// time, in milliseconds, between state transitions
#define TL_FLL_SWITCH_TIME 500
#define TL_PLL_TIME 1000
// must be greater than TL_PLL_TIME
// in the interval between TL_PLL_TIME and TL_PULLIN_TIME the data for detecting the nav bit edge is gathered
#define TL_PULLIN_TIME 1500

#define DOPPLER_RADIUS 5000
#define NUM_COARSE_DOPPLERS 101  // add one for center, keep sides even
#define FINE_DOPPLER_RADIUS (DOPPLER_RADIUS * 2 / NUM_COARSE_DOPPLERS)
#define NUM_FINE_DOPPLERS 251

#define CODE_FREQ      1.023e6
#define CHIPS_PER_CODE 1023

#define CARRIER_FREQ   1.57542e9
#define CARRIER_AID_SF (CODE_FREQ / CARRIER_FREQ)
#define CORRELATOR_BUF_SIZE (2000)


//#define SAMPLING_FREQ (double)38.192e6 // for bladeRF
//#define SAMPLING_FREQ (double)40.0e6
#define SAMPLING_FREQ (double)53.0e6 // for NTlab

#define ts 1/SAMPLING_FREQ;
#define CODE_TIME_INC (1/SAMPLING_FREQ)



//#define BUFFLENGTH (ACQ_MS*(SAMPLING_FREQ/1000)) //samplesPerCode
//#define BUFFLENGTH 1600000
//#define BUFFLENGTH 440000
#define BUFFLENGTHint (int)(ACQ_MS*(SAMPLING_FREQ/1000))
#define BUFFLENGTH (int)(11*(SAMPLING_FREQ/1000))
//#define BUFFLENGTH 53000
#define SAMPLES_MS (SAMPLING_FREQ/1000)

#define IF (double)10e6// for bladeRF

#define doppler_max = 5000
#define doppler_min = -5000
#define doppler_step = 100
//#define doppler_bins = (doppler_max-doppler_min)/doppler_step
#define doppler_bins = 100

//navigation data
#define MS_PER_NAV_BIT 20
#define SUBFRAME_LENGTH 300
#define PREAMBLE_LENGTH 8
#define MAX_PREAMBLE_CANDIDATES 10
#define PAYLOAD_LENGTH 30
#define NAV_DATA_BIT_DURATION 0.02


#define TWO_NEG_5				0.03125
#define TWO_NEG_29			1.86264514923095703125e-9
#define TWO_NEG_31			4.656612873077392578125e-10
#define TWO_NEG_33			1.16415321826934814453125e-10
#define TWO_NEG_43			1.136868377216160297393798828125e-13
#define TWO_NEG_19			1.9073486328125e-6
#define TWO_NEG_55			2.77555756e-17
#define TWO_POS_4				16

//Pseudo
#define NAV_C						299792458
#define NAV_OMEGAE_DOT	7.2921151467e-005
#define NAV_C						299792458
#define NAV_GM					3.986005e14
#define NAV_F						-4.442807633e-010
#define NAV_DATA_BIT_DURATION	0.02
#define R2D 57.295779513082320876798154814105           /* radians to degrees */
#define D2R 1.0/R2D                                    /* degrees to radians */

//Processing
#define F          ( -4.44280763310e-10 ) //constant to calculate relativistic correction term
#define C_         ( 2.99792458e8 ) //speed of lite
#define mu         ( 3.986005e14) //universal gravitational parameter
#define sqrt_of_mu ( 1.996498184322e7) //square root from mu value
#define w_earth    ( 7.2921151467e-5 )//Earth's rotation rate
#define R_e        ( 6371302 )//average radius of the Earth;
#define Pi         ( M_PI )//pi value
#define a_WGS      ( 6378137 )//semi-major axis of WGS-84 ellipsoid
#define f_WGS      ( 1/298.257223563 )//flattening factor of WGS-84 ellipsoid
#define f_L1       (1575.42e6)//GPS frequency of L1 band
#define f_L2       (1227.60e6)//GPS frequency of L2 band

#define GPSf1s								(f_L1*f_L1)
#define GPSf2s								(f_L2*f_L2)
#define GPSf12s								(GPSf1s-GPSf2s)
#define GPSf1ION							(GPSf2s/GPSf12s)
#define GPSf2ION							(GPSf1s/GPSf12s)

