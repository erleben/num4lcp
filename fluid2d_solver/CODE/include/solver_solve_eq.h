#ifndef SOLVER_SOLVE_EQ_H
#define SOLVER_SOLVE_EQ_H

#include <util_params.h>

#include <gpu_print_versions.h>
#include <gpu_print_current_device.h>
#include <gpu_check_error.h>
#include <gpu_safe_call.h>
#include <gpu_cpu_timer.h>
#include <gpu_cuda_timer.h>

#include <solver_profiling_data.h>
#include <solver_linear_system_solver.h>
#include <solver_status_codes.h>
#include <solver_make_coo_matrix.h>

#include <cusp/coo_matrix.h>
#include <cusp/csr_matrix.h>
#include <cusp/linear_operator.h>
#include <cusp/precond/ainv.h>
#include <cusp/precond/diagonal.h>
#include <cusp/monitor.h>
#include <cusp/io/matrix_market.h>
#include <cusp/convert.h>

#include <iosfwd>
#include <cassert>
#include <vector>


/**
 *
 * @tparam M2    The target memory space can be either cusp::host_memory or
 *               cusp::device_memory. This essesially controls whether a GPU
 *               or CPU computation is done.
 */
template<typename M2>
inline __host__ void solve_eq(
                       unsigned int n
                       , unsigned int k
                       , std::vector<unsigned int> const & A_row_indices
                       , std::vector<unsigned int> const & A_column_indices
                       , std::vector<float       > const & A_values
                       , std::vector<float       >       & x_values
                       , std::vector<float       > const & b_values
                       , Params                    const & params
                       , ProfilingData                   & profiling_data
                       )
{
  typedef int                 I;
  typedef float               T;
  typedef cusp::host_memory   M1;

  CHECK_ERROR("solve_eq invoked");

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
    assert( num_devices > device_number || !"solve_eq(): Invalid device number: " + device_number);

    if (num_devices-1 < device_number)
    {
      std::cout << "solve_eq(): Invalid device ID, device ID was "
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

    assert( byte_needed < props.totalGlobalMem || !"solve_eq(): Not enough GPU memory");
  }

  //--- Create host COO matrices and vectors
  cusp::coo_matrix<I,  T,  M1>  A1 = make_coo_matrix<I, T, M1>( n
                                                              , k
                                                              , A_row_indices
                                                              , A_column_indices
                                                              , A_values
                                                              );
  cusp::array1d<T, M1> x1(x_values);
  cusp::array1d<T, M1> b1(b_values);
  CHECK_ERROR("Host matrix vectors created and initialized");

  //--- Set stopping criteria
  unsigned int desired_max_iterations = params.eq_max_iterations();
  unsigned int necessary_max_iterations = n / 10;

  unsigned int const eq_max_iterations      = max(desired_max_iterations, necessary_max_iterations );
  T            const eq_absolute_tolerance  = params.eq_relative_tolerance();
  T            const eq_relative_tolerance  = params.eq_relative_tolerance();
  unsigned int const eq_choice              = params.sub_solver_choice();

  cusp::monitor<T> monitor(
                           b1
                           , eq_max_iterations
                           , eq_absolute_tolerance
                           , eq_relative_tolerance
                           );

  CHECK_ERROR("Monitor created");


  if(params.profiling() && params.record_time())
  {
    cpu_timer.stop();
    profiling_data.m_host_initialization_time = cpu_timer();
    gpu_timer.start();
  }
  //--- Convert COO‐>CSR on the host and transfer to the device ----------------
  cusp::csr_matrix<I, T,  M2> A2 = A1;
  cusp::array1d<T, M2>        x2 = x1;
  cusp::array1d<T, M2>        b2 = b1;
  CHECK_ERROR("Device matrix vectors created and initialized");

  if(params.profiling() && params.record_time())
  {
    gpu_timer.stop();
    profiling_data.m_host_to_device_time = gpu_timer();
    gpu_timer.start();
  }

  switch(params.preconditioner_choice())
  {
    case 0:
    {
      cusp::identity_operator<T,M2, I> P2(n, n);
      
      linear_system_solver(
                           A2
                           , b2
                           , x2
                           , monitor
                           , P2
                           , params.sub_solver_choice()
                           );
    }
      break;
    case 1:
    {
      cusp::precond::diagonal<T,M2> P2(A2);

      linear_system_solver(
                           A2
                           , b2
                           , x2
                           , monitor
                           , P2
                           , params.sub_solver_choice()
                           );
    }
      break;
    case 2:
    {
      cusp::precond::scaled_bridson_ainv<T,M2> P2(A2, .1);

      linear_system_solver(
                           A2
                           , b2
                           , x2
                           , monitor
                           , P2
                           , params.sub_solver_choice()
                           );
    }
      break;
    case 3:
    {
      cusp::precond::scaled_bridson_ainv<T,M2> P2(A2, 0, 10);

      linear_system_solver(
                           A2
                           , b2
                           , x2
                           , monitor
                           , P2
                           , params.sub_solver_choice()
                           );
    }
      break;
    case 4:
    {
      cusp::precond::bridson_ainv<T,M2> P2(A2, 0, -1, true, 2);

      linear_system_solver(
                           A2
                           , b2
                           , x2
                           , monitor
                           , P2
                           , params.sub_solver_choice()
                           );
    }
      break;
  };
  CHECK_ERROR("Solver finished");

  if(params.profiling() && params.record_time())
  {
    gpu_timer.stop();
    profiling_data.m_computation_time = gpu_timer();
  }

  if(params.profiling() && params.record_time())
  {
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

  if(params.profiling() && params.record_convergence())
  {
    profiling_data.m_eq_iteration = monitor.iteration_count();
    cusp::copy(monitor.residuals, profiling_data.m_eq_residuals);
    profiling_data.m_eq_status    = monitor.converged() ? EQ_SOLVER_OK : EQ_SOLVER_FAILED;
  }
}


__host__ void solve_eq_host(
                            unsigned int n
                            , unsigned int k
                            , std::vector<unsigned int> const & A_row_indices
                            , std::vector<unsigned int> const & A_column_indices
                            , std::vector<float       > const & A_values
                            , std::vector<float       >       & x_values
                            , std::vector<float       > const & b_values
                            , Params                    const & params
                            , ProfilingData                   & profiling_data
                            )
{
  typedef cusp::host_memory   M2;
  
  solve_eq<M2>(n
               , k
               , A_row_indices
               , A_column_indices
               , A_values
               , x_values
               , b_values
               , params
               , profiling_data
               );
}


__host__ void solve_eq_device(
                              unsigned int n
                              , unsigned int k
                              , std::vector<unsigned int> const & A_row_indices
                              , std::vector<unsigned int> const & A_column_indices
                              , std::vector<float       > const & A_values
                              , std::vector<float       >       & x_values
                              , std::vector<float       > const & b_values
                              , Params                    const & params
                              , ProfilingData                   & profiling_data
                              )
{
  typedef cusp::device_memory M2;
  
  solve_eq<M2>(n
               , k
               , A_row_indices
               , A_column_indices
               , A_values
               , x_values
               , b_values
               , params
               , profiling_data
               );
}


// SOLVER_SOLVE_EQ_H
#endif
