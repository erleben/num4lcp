#ifndef GPU_GET_DEVICE_COUNT_H
#define GPU_GET_DEVICE_COUNT_H

#include <gpu_safe_call.h>

#include <string>

/**
 * Get the device count.
 */
inline __host__ int get_device_count()
{
  int count = 0;
  SAFE_CALL( cudaGetDeviceCount( &count) );
  return count;
}

// GPU_GET_DEVICE_COUNT_H
#endif
