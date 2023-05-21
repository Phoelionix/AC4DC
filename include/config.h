/**
 * @file config.h
 * @brief 
 */

#pragma once

const int GLOBAL_BSPLINE_ORDER = 3;  // 1 = rectangles, 2=linear, 3=quadratic  ATTENTION I have not yet gone through and ensured that there are no local redefinitions of the order.

#ifdef DEBUG
// #define OUTPUT_TBR_TO_CERR
// #define OUTPUT_DQDT_TO_CERR
// #define OUTPUT_DFDT_TO_CERR
#endif

//  #define NO_TBR
//  #define NO_EE #TODO try turning off again
//  #define NO_EII

//#define DEBUG_BOUND
//#define SWITCH_OFF_DYNAMIC_BOUNDS
//#define SWITCH_OFF_ALL_DYNAMIC_UPDATES  // should be equivalent (or almost equivalent, I can't remember if gaussian quadrature will ruin us here) to switching off dynamic bounds if working properly