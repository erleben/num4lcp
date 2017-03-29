#ifndef SOLVER_SOLVE_LCP_H
#define SOLVER_SOLVE_LCP_H

#include <gpu_print_versions.h>
#include <gpu_print_current_device.h>
#include <gpu_check_error.h>
#include <gpu_safe_call.h>
#include <gpu_cpu_timer.h>
#include <gpu_cuda_timer.h>

#include <util_params.h>

#include <solver_min_map_newton.h>
#include <solver_psor.h>
#include "solver_gradient_descent.h"
#include <solver_profiling_data.h>
#include <solver_make_coo_matrix.h>

#include <cusp/coo_matrix.h>
#include <cusp/csr_matrix.h>
#include <cusp/linear_operator.h>
#include <cusp/precond/ainv.h>
#include <cusp/precond/diagonal.h>

#include <thrust/transform.h>
#include <thrust/functional.h>

#include <iosfwd>
#include <cassert>
#include <vector>


template<typename M2>
inline __host__ void solve_lcp(
                                      unsigned int n
                                      , unsigned int k
                                      , std::vector<unsigned int> const & A_row_indices
                                      , std::vector<unsigned int> const & A_column_indices
                                      , std::vector<float       > const & A_values
                                      , std::vector<float       >       & x_values
                                      , std::vector<float       > const & b_values
                                      , std::vector<bool        > const & is_bilateral
                                      , Params                    const & params
                                      , ProfilingData                   & profiling_data
                                      )
{
  typedef int                  I;
  typedef float                T;
  typedef cusp::host_memory    M1;

  CHECK_ERROR("solve_lcp_device invoked");
  
  CPUTimer  cpu_timer;
  CUDATimer gpu_timer;

  if(params.profiling() && params.record_time())
  {
    cpu_timer.start();
  }

  if(params.use_gpu())
  {
    int device_number = params.device_number();

    //--- Check number of available GPU devices
    //--- and do a memory check

    int num_devices = 0;
    
    SAFE_CALL( cudaGetDeviceCount(&num_devices) );

    // Device numbering is zero indexed, so device 1 has id 0. So we have to
    // subtract one from number of devices
    assert( num_devices > device_number || !"solve_lcp_device(): Invalid device number: " + device_number);
    
    if (num_devices-1 < device_number)
    {
      std::cout << "solve_lcp_device(): Invalid device ID, device ID was "
                << device_number
                << " of ["
                << num_devices-1
                << "] available."
                << std::endl;
      
      exit(EXIT_FAILURE);
    }
    
    cudaDeviceProp props;
    
    SAFE_CALL( cudaGetDevice( &device_number) );
    SAFE_CALL( cudaGetDeviceProperties(&props, device_number) );
    
    unsigned int const byte_needed = n*(5*6 + 2);
    
    assert( byte_needed < props.totalGlobalMem || !"solve_lcp_device(): Not enough GPU memory");
  }
  
  //--- Allocate space and create host matrices and vectors
  cusp::coo_matrix<I,  T,  M1>  A1coo = make_coo_matrix<I, T, M1>( n
                                                                 , k
                                                                 , A_row_indices
                                                                 , A_column_indices
                                                                 , A_values
                                                                 );

  cusp::csr_matrix<I, T,  M1>  A1 = A1coo;
  cusp::array1d<T, M1> x1(x_values);
  cusp::array1d<T, M1> b1(b_values);
  cusp::array1d<bool, M1> is_bilateral1(is_bilateral);
  CHECK_ERROR("Host matrices and vectors created and initialized");

  //--- LCP is defined as y = A x + b, but we are feeding in A x = b so we
  //--- negate right hand side vector.
  thrust::transform( b_values.begin()
                    , b_values.end()
                    , b1.begin()
                    , thrust::negate<T>()
                    );

  if(params.profiling() && params.record_time())
  {
    cpu_timer.stop();
    profiling_data.m_host_initialization_time = cpu_timer();
    gpu_timer.start();
  }

  if(params.psor_warmstart())
  {
    if(params.profiling() && params.record_time())
    {
      cpu_timer.start();
    }

    psor(
         A1
         , b1
         , x1
         , is_bilateral1
         , params
         , profiling_data.m_psor_status
         , profiling_data.m_psor_iteration
         , profiling_data.m_psor_residuals
         );

    if(params.profiling() && params.record_time())
    {
      cpu_timer.stop();
      profiling_data.m_psor_time = cpu_timer();
    }
  }
  else if (params.grad_warmstart())
  {
    if(params.profiling() && params.record_time())
    {
      cpu_timer.start();
    }
      gradient_descent(
              A1, b1, x1, is_bilateral1, params, profiling_data.m_grad_status, profiling_data.m_grad_iteration,
              profiling_data.m_grad_residuals
      );

    if(params.profiling() && params.record_time())
    {
      cpu_timer.stop();
    }
  }

  if(params.profiling() && params.record_time())
  {
    gpu_timer.start();
  }

  //--- Transfer matrix and vecotors to the device
  cusp::csr_matrix<I, T,  M2> A2            = A1;
  cusp::array1d<T, M2>        x2            = x1;
  cusp::array1d<T, M2>        b2            = b1;
  cusp::array1d<bool, M2>     is_bilateral2 = is_bilateral1;
  CHECK_ERROR("Device matrix and vectors created");

  if(params.profiling() && params.record_time())
  {
    gpu_timer.stop();
    profiling_data.m_host_to_device_time = gpu_timer();
    gpu_timer.start();
  }
  
  //--- Solve the LCP ----------------------------------------------------------
  switch(params.preconditioner_choice())
  {
    case 0:
    {
      cusp::identity_operator<T,M2,I> P2(n, n);

      min_map_newton(
                     A2
                     , b2
                     , x2
                     , is_bilateral2
                     , P2
                     , params
                     , profiling_data.m_newton_status
                     , profiling_data.m_newton_iteration
                     , profiling_data.m_newton_residuals
                     );
    }
      break;
    case 1:
    {
      cusp::precond::diagonal<T,M2> P2(A2);

      min_map_newton(
                     A2
                     , b2
                     , x2
                     , is_bilateral2
                     , P2
                     , params
                     , profiling_data.m_newton_status
                     , profiling_data.m_newton_iteration
                     , profiling_data.m_newton_residuals
                     );
    }
      break;
    case 2:
    {
      cusp::precond::scaled_bridson_ainv<T,M2> P2(A2, .1);

      min_map_newton(
                     A2
                     , b2
                     , x2
                     , is_bilateral2
                     , P2
                     , params
                     , profiling_data.m_newton_status
                     , profiling_data.m_newton_iteration
                     , profiling_data.m_newton_residuals
                     );
    }
      break;
    case 3:
    {
      cusp::precond::scaled_bridson_ainv<T,M2> P2(A2, 0, 10);

      min_map_newton(
                     A2
                     , b2
                     , x2
                     , is_bilateral2
                     , P2
                     , params
                     , profiling_data.m_newton_status
                     , profiling_data.m_newton_iteration
                     , profiling_data.m_newton_residuals
                     );
    }
      break;
    case 4:
    {
      cusp::precond::bridson_ainv<T,M2> P2(A2, 0, -1, true, 2);

      min_map_newton(
                     A2
                     , b2
                     , x2
                     , is_bilateral2
                     , P2
                     , params
                     , profiling_data.m_newton_status
                     , profiling_data.m_newton_iteration
                     , profiling_data.m_newton_residuals
                     );
    }
      break;
  };
  CHECK_ERROR("Solver finished");
  
  if(params.profiling() && params.record_time())
  {
    gpu_timer.stop();
    profiling_data.m_computation_time = gpu_timer();
    gpu_timer.start();
  }
  
  //--- Read back solution from device to host ---------------------------------
  x1 = x2;
  CHECK_ERROR("Readback solution from device to host");
  
  thrust::copy( x1.begin(), x1.end(), x_values.begin() );
  
  if(params.profiling() && params.record_time())
  {
    gpu_timer.stop();
    profiling_data.m_device_to_host_time = gpu_timer();
  }
}


inline __host__ void solve_lcp_host(
                               unsigned int n
                               , unsigned int k
                               , std::vector<unsigned int> const & A_row_indices
                               , std::vector<unsigned int> const & A_column_indices
                               , std::vector<float       > const & A_values
                               , std::vector<float       >       & x_values
                               , std::vector<float       > const & b_values
                               , std::vector<bool        > const & is_bilateral
                               , Params                    const & params
                               , ProfilingData                   & profiling_data
                               )
{
  typedef cusp::host_memory M2;

  solve_lcp<M2>(
                n
                , k
                , A_row_indices
                , A_column_indices
                , A_values
                , x_values
                , b_values
                , is_bilateral
                , params
                , profiling_data
                );
}


inline __host__ void solve_lcp_device(
                               unsigned int n
                               , unsigned int k
                               , std::vector<unsigned int> const & A_row_indices
                               , std::vector<unsigned int> const & A_column_indices
                               , std::vector<float       > const & A_values
                               , std::vector<float       >       & x_values
                               , std::vector<float       > const & b_values
                               , std::vector<bool        > const & is_bilateral
                               , Params                    const & params
                               , ProfilingData                   & profiling_data
                               )
{
  typedef cusp::device_memory M2;

  solve_lcp<M2>(
                n
                , k
                , A_row_indices
                , A_column_indices
                , A_values
                , x_values
                , b_values
                , is_bilateral
                , params
                , profiling_data
                );
}


// SOLVER_SOLVE_LCP_H
#endif
