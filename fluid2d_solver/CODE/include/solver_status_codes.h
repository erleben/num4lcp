#ifndef SOLVER_STATUS_CODES_H
#define SOLVER_STATUS_CODES_H

/**
 * Error/status return codes.
 */

unsigned int const ITERATING             = 0u;
unsigned int const NON_DESCEND_DIRECTION = 1u;
unsigned int const STAGNATION            = 2u;
unsigned int const RELATIVE_CONVERGENCE  = 3u;
unsigned int const ABSOLUTE_CONVERGENCE  = 4u;
unsigned int const LOCAL_MINIMA          = 5u;
unsigned int const MAX_LIMIT             = 6u;
unsigned int const LINE_SEARCH_FAILURE   = 7u;
unsigned int const SUFFICIENT_DECREASE   = 8u;
unsigned int const TOO_SMALL_STEP        = 9u;
unsigned int const SUB_SOLVER_FAILURE    = 10u;
unsigned int const EQ_SOLVER_OK          = 11u;
unsigned int const EQ_SOLVER_FAILED      = 12u;


// SOLVER_STATUS_CODES_H
#endif
