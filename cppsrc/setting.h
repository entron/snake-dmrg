#ifndef SETTING_H
#define SETTING_H


#define NUMBER_OF_GOOD_QUANTUM_NUMBER 1
#define NUMBER_OF_KINDS_OF_PARTICLES 1
#define FERMIONSIGN 0

///Whether to calculate all sites.1 means yes, and 0 means no.
#define CALCULATE_ALL_SITES 0

#define USE_INFINITE_DMRG 1


#define USE_EVOLVE 1
///The following four options are only valid when USE_EVOLVE
  #define TARGET_TWO_WF 0 
  #define CAL_DURING_TEMPERATURE_EVOLVE 1
  #define CAL_DURING_TIME_EVOLVE 1



#define   Max_Truncate_Error 1e-12


#include "gqnmat.h"
#define LA_COMPLEX_SUPPORT
#include <gmc.h>
#include <gmd.h>
typedef snake::math::GQNMat<LaGenMatDouble>  Rmatrix;
typedef snake::math::GQNMat<LaGenMatComplex>   Cmatrix;


#endif
