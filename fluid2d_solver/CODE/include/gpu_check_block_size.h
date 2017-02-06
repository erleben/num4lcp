#ifndef GPU_CHECK_BLOCK_SIZE_H
#define GPU_CHECK_BLOCK_SIZE_H

#include <gpu_safe_call.h>

#include <iostream>

/**
 * Based on the GPU specs and the register count for a kernel, the best block size is selected and
 * compared with the current block size.
 *
 * @param func     A (__global__) functor
 * @param int      The current block size
 * @return         The best block size.
 */
template<typename T>
inline __host__ int check_block_size( T const & func, int block_size)
{
  cudaFuncAttributes attrib;
  cudaDeviceProp props;
  
  // Get current device
  int device_number = 0;
  SAFE_CALL( cudaGetDevice( &device_number) );
  
  // The largest block size is preferred, since maximum 8 blocks can be active
  // in a single SM.
  size_t best_block_size = block_size;
  
  // Retrieve info
  SAFE_CALL( cudaFuncGetAttributes( &attrib, func ) );     
  SAFE_CALL( cudaGetDeviceProperties(&props, device_number) );
  
  int blocks;
  int wasted = attrib.maxThreadsPerBlock;
  for (size_t test_bs = 32u; test_bs < attrib.maxThreadsPerBlock; test_bs += 32u) 
  {
    blocks = attrib.maxThreadsPerBlock / test_bs;
    // block count must be 8 or below.
    if (blocks <= 8) 
    {
      // Reduce wasted threads?
      if( (attrib.maxThreadsPerBlock % test_bs) < wasted) 
      { 
        best_block_size = test_bs;
        wasted = attrib.maxThreadsPerBlock % test_bs;
      }
    }
  }
  
  blocks = (attrib.maxThreadsPerBlock / block_size);
  int warps = (blocks*block_size)/32;
  
  std::cout << "Best block size    :" << best_block_size                             << std::endl;
  std::cout << "Current block size :" << block_size                                  << std::endl;
  std::cout << "Wasted threads     :" << attrib.maxThreadsPerBlock % block_size      << std::endl;
  std::cout << "Active blocks      :" << blocks                                      << std::endl;
  std::cout << "Active warps       :" << warps                                       << std::endl;
  std::cout << "Intermidiate clock cycles required to hide latency:" << 200/warps    << std::endl;
  std::cout << std::endl;
  
  return best_block_size;
}

// GPU_CHECK_BLOCK_SIZE_H
#endif
