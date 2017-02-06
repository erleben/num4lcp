#ifndef GPU_PRINT_VERSIONS_H
#define GPU_PRINT_VERSIONS_H

#include <cuda.h>

#include <thrust/version.h>
#include <cusp/version.h>

#include <sstream>
#include <string>

inline __host__ std::string print_versions()
{
  std::stringstream output;

  int const cuda_major =  CUDA_VERSION / 1000;
  int const cuda_minor = (CUDA_VERSION % 1000) / 10;
  
  output << "CUDA   v" << cuda_major   << "." << cuda_minor   << std::endl;
  
  int const thrust_major = THRUST_MAJOR_VERSION;
  int const thrust_minor = THRUST_MINOR_VERSION;
  int const thrust_subminor = THRUST_SUBMINOR_VERSION;

  output << "Thrust v" << thrust_major << "." << thrust_minor << "." << thrust_subminor << std::endl;
  
  int const cusp_major    = CUSP_MAJOR_VERSION;
  int const cusp_minor    = CUSP_MINOR_VERSION;
  int const cusp_subminor = CUSP_SUBMINOR_VERSION;
  
  output << "Cusp   v" << cusp_major << "." << cusp_minor << "." << cusp_subminor << std::endl;

  output.flush();

  return output.str();
}


// GPU_PRINT_VERSIONS_H
#endif
