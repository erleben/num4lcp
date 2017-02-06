#ifndef GPU_PRINT_FUNCTION_H
#define GPU_PRINT_FUNCTION_H

#include <gpu_safe_call.h>

#include <iostream>

/**
 * Print stats of a GPU function
 * @param func     A (__global__) functor
 * @param name     A human readable name to identify output by,
 */
template<typename T>
inline __host__ void print_function(  T const & func , char const * name)
{
  cudaFuncAttributes attrib;
  
  SAFE_CALL( cudaFuncGetAttributes( &attrib, func ) ); 
  
  std::cout << "Function                           :" << name << std::endl;
  std::cout << "Size of constant memory in bytes   :" << attrib.constSizeBytes     << std::endl;
  std::cout << "Size of local memory in bytes      :" << attrib.localSizeBytes     << std::endl;
  std::cout << "Maximum number of threads per block:" << attrib.maxThreadsPerBlock << std::endl;
  std::cout << "Number of registers used           :" << attrib.numRegs            << std::endl;
  std::cout << "Size of shared memory in bytes     :" << attrib.sharedSizeBytes    << std::endl;
  std::cout << std::endl;
}


// GPU_PRINT_FUNCTION_H
#endif
