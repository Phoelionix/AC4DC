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

//#define NO_TBR
//#define NO_EE   // This can break the dynamic grid late in the simulation at low energies.
//#define NO_EII

//#define NO_MINISTEPS
#define NO_MINISTEP_UPDATING  // ministep updates not working atm, possibly as lagrange polynomial is a bad thing to use here and I should just use linear interpolation.

//#define DEBUG_BOUND
//#define SWITCH_OFF_DYNAMIC_BOUNDS
//#define SWITCH_OFF_ALL_DYNAMIC_UPDATES  // should be equivalent (or almost equivalent, I can't remember if gaussian quadrature will make a difference here) to switching off dynamic bounds if working properly

//#define NO_PLOTTING
//#define NO_BACKUP_SAVING
/////////////


#ifdef NO_MINISTEPS 
    #ifndef NO_MINISTEP_UPDATING
    #define NO_MINISTEP_UPDATING
    #endif
#endif