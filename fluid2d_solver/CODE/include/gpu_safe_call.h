#ifndef GPU_SAFE_CALL_H
#define GPU_SAFE_CALL_H

#include <cuda_runtime.h>    // needed for CUDA C++ runtime api
#include <cstdlib>           // needed for exit() function
#include <iostream>          // needed for std::cerr

/**
 * @file  This file contains a small collection of macros that makes it easier to work with the CUDA runtime library.
 * The functionality has been greatly inspired by the CUDA Utility (CUTIL) library that is used in the CUDA SDK.
 * However, since CUTIL is not part of the CUDA SDK we have for reason of self-containness borrowed a few ideas
 * from CUTIL and reimplemented them here in our own "flavour".  Should CUTIL ever become part of the CUDA SDK
 * then it makes sense to replace our functionality with equivalents CUTIL counter parts.
 */

/**
 * Cuda Runtime Safe Call Macro.
 * This macro automaticall tests the return code from a cuda runtime call and generates
 * an error message if the call is not a success and exits program execution.
 *
 * @param call       The cuda runtime call. 
 */
#  define SAFE_CALL( call) {                                                 \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        std::cerr << "Cuda error in file "                                   \
                  << __FILE__                                                \
                  <<  " line "                                               \
                  << __LINE__                                                \
                  << " : "                                                   \
                  <<  cudaGetErrorString( err)                               \
                  << std::endl;                                              \
        exit(EXIT_FAILURE);                                                  \
} }

// GPU_SAFE_CALL_H
#endif 
