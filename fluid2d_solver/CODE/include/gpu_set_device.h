#ifndef GPU_SET_DEVICE_H
#define GPU_SET_DEVICE_H

#include <gpu_safe_call.h>

/**
 * Set the device.
 */
inline __host__ void set_device(int const & device_number)
{
  SAFE_CALL( cudaSetDevice( device_number ) );
}

// GPU_SET_DEVICE_H
#endif
