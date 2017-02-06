#ifndef GPU_PRINT_DEVICE_H
#define GPU_PRINT_DEVICE_H

#include <gpu_safe_call.h>

#include <sstream>
#include <string>

/**
 * Print device properties onto screen.
 *
 * @param device_number    The number of the device whos properties are printed.
 */
inline __host__ std::string print_device(int device_number)
{
  std::stringstream output;

  cudaDeviceProp props;
  SAFE_CALL( cudaGetDeviceProperties(&props, device_number) );
  output << "Device "<< device_number << ":" <<  props.name              << std::endl;
  output << "Global memory      :" <<  props.totalGlobalMem              << std::endl;
  output << "Shared memory      :" <<  props.sharedMemPerBlock           << std::endl;
  output << "Registers          :" <<  props.regsPerBlock                << std::endl;
  output << "Warp               :" <<  props.warpSize                    << std::endl;
  output << "Pitch              :" <<  props.memPitch                    << std::endl;
  output << "Threads per block  :" <<  props.maxThreadsPerBlock          << std::endl;
  output << "Clock frequency    :" <<  props.clockRate                   << std::endl;
  output << "Constant memory    :" <<  props.totalConstMem               << std::endl;
  output << "Compute capability :" <<  props.major << "."  <<  props.minor  << std::endl;
  output << "Texture alignment  :" <<  props.textureAlignment            << std::endl;
  output << "Concurrently copy  :" << (props.deviceOverlap == 1 )        << std::endl;
  output << "Multiprocessors    :" <<  props.multiProcessorCount         << std::endl;
  output << "Run time limit     :" << (props.kernelExecTimeoutEnabled==1) << std::endl;
  output << "Integrated GPU     :" << (props.integrated==1)              << std::endl;
  output << "Map host memory    :" << (props.canMapHostMemory==1)        << std::endl;
  output << "Compute_mode       :";
  switch( props.computeMode)
  {
    case cudaComputeModeDefault:
      output << "Default mode" << std::endl;
      break;
    case cudaComputeModeExclusive: 
      output << "Compute-exclusive mode"<< std::endl;
      break;
    case cudaComputeModeProhibited: 
      output << "Compute-prohibited mode"<< std::endl;
      break;
  };     
  output << "Block Dim          :" <<  props.maxThreadsDim[0] << " X " << props.maxThreadsDim[1] << " X " << props.maxThreadsDim[2] << std::endl;
  output << "Grid Dim           :" <<  props.maxGridSize[0]   << " X " << props.maxGridSize[1]   << " X " << props.maxGridSize[2]   << std::endl;
  output << std::endl;

  output.flush();

  return output.str();
}

// GPU_PRINT_DEVICE_H
#endif
