#ifndef GPU_SHARED_MEMORY_H
#define GPU_SHARED_MEMORY_H

#include <cuda_runtime.h>    // needed for CUDA C++ runtime api

/**
 * Utility class used to avoid linker errors with extern
 * unsized shared memory arrays with templated type
 *
 * Example:
 *
 * T * myarray = SharedMemory<T>();
 */
template<class T>
struct SharedMemory
{
  __device__ inline operator       T*()
  {
    extern __shared__ int __smem[];
    return (T*)__smem;
  }
  
  __device__ inline operator const T*() const
  {
    extern __shared__ int __smem[];
    return (T*)__smem;
  }
};

/**
 * specialize for float to avoid unaligned memory 
 * access compile errors
 */
template<>
struct SharedMemory<float>
{
  __device__ inline operator       float*()
  {
    extern __shared__ float __smem_d[];
    return (float*)__smem_d;
  }
  
  __device__ inline operator const float*() const
  {
    extern __shared__ float __smem_d[];
    return (float*)__smem_d;
  }
};


// GPU_SHARED_MEMORY_H
#endif 
