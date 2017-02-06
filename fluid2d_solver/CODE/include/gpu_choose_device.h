#ifndef GPU_CHOOSE_DEVICE_H
#define GPU_CHOOSE_DEVICE_H

#include <gpu_safe_call.h>

/**
 * Scans through all GPU devices on the system and picks the device that is best suited to do the GPU computations on.
 *
 * @return    The device number.
 */
inline __host__ int choose_device()
{
  int best_gpu       = 0;
  int number_of_devices = 0;
  SAFE_CALL( cudaGetDeviceCount(&number_of_devices) );
  if (number_of_devices > 1)
  {
    int max_processors = 0;
    for (int device_number = 0; device_number < number_of_devices; ++device_number)
    {
      cudaDeviceProp props;
      SAFE_CALL( cudaGetDeviceProperties(&props, device_number) );
      if (props.multiProcessorCount > max_processors)
      {
        best_gpu = device_number;
        max_processors = props.multiProcessorCount; 
      }
    }
  } 
  return best_gpu;
}

// GPU_CHOOSE_DEVICE_H
#endif
