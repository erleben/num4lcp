#ifndef GPU_CHECK_ERROR_H
#define GPU_CHECK_ERROR_H

#include <cuda_runtime.h>    // needed for CUDA C++ runtime api
#include <cstdlib>           // needed for exit() function
#include <iostream>          // needed for std::cerr

/**
 * Cuda Runtime Check for Errors Macro.
 * This macro is usefull if one wants to verify that no errors have occured in the cuda
 * runtime. One can specify a text message that will be printed in case of an error. The
 * macro will exit program if an error is detected.
 *
 * @param errorMessage    The text to be printed in case of an error.
 */
#ifdef _DEBUG
#  define CHECK_ERROR(errorMessage) {                                        \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        std::cerr << "Cuda error in file "                                   \
                  << __FILE__                                                \
                  <<  " line "                                               \
                  << __LINE__                                                \
                  << " : "                                                   \
                  <<  cudaGetErrorString( err)                               \
                  << std::endl;                                              \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    err = cudaThreadSynchronize();                                           \
    if( cudaSuccess != err) {                                                \
        std::cerr << "Cuda error in file "                                   \
                  << __FILE__                                                \
                  <<  " line "                                               \
                  << __LINE__                                                \
                  << " : "                                                   \
                  <<  cudaGetErrorString( err)                               \
                  << std::endl;                                              \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
}
#else
#  define CHECK_ERROR(errorMessage) {                                        \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        std::cerr << "Cuda error in file "                                   \
                  << __FILE__                                                \
                  <<  " line "                                               \
                  << __LINE__                                                \
                  << " : "                                                   \
                  <<  cudaGetErrorString( err)                               \
                  << std::endl;                                              \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
}
#endif


// GPU_CHECK_ERROR_H
#endif 
