#ifndef GPU_PRINT_CURRENT_DEVICE_H
#define GPU_PRINT_CURRENT_DEVICE_H

#include <gpu_safe_call.h>
#include <gpu_print_device.h>

#include <string>

/**
 * Print current device properties onto screen.
 */  
inline __host__ std::string print_current_device()
{
  int device_number = 0;
  SAFE_CALL( cudaGetDevice( &device_number) );
  return print_device( device_number );
}

// GPU_PRINT_CURRENT_DEVICE_H
#endif
