#ifndef SOLVER_PROFILING_DATA_H
#define SOLVER_PROFILING_DATA_H

#include <vector>
#include <cusp/array1d.h>

template<typename T>
inline std::vector<T> make_std_vector(cusp::array1d<T,cusp::host_memory> const & x)
{
  std::vector<T> y;

  y.resize(x.size() );

  thrust::copy( x.begin(), x.end(), y.begin() );

  return y;
}


class ProfilingData
{
public:

  float m_host_initialization_time;
  float m_host_to_device_time;
  float m_device_to_host_time;
  float m_computation_time;
  float m_psor_time;

  unsigned int m_newton_status;
  unsigned int m_psor_status;
  unsigned int m_eq_status;
  unsigned int m_grad_status;

  unsigned int m_newton_iteration;
  unsigned int m_psor_iteration;
  unsigned int m_eq_iteration;
  unsigned int m_grad_iteration;

  cusp::array1d<float, cusp::host_memory> m_newton_residuals;
  cusp::array1d<float, cusp::host_memory> m_psor_residuals;
  cusp::array1d<float, cusp::host_memory> m_eq_residuals;
  cusp::array1d<float, cusp::host_memory> m_grad_residuals;

public:

  ProfilingData()
    : m_host_initialization_time(0.0)
    , m_host_to_device_time(0.0)
    , m_device_to_host_time(0.0)
    , m_computation_time(0.0)
    , m_psor_time(0.0)
    , m_newton_status(0u)
    , m_psor_status(0u)
    , m_eq_status(0u)
    , m_grad_status(0u)
    , m_newton_iteration(0u)
    , m_psor_iteration(0u)
    , m_eq_iteration(0u)
    , m_grad_iteration(0u)
    , m_newton_residuals()
    , m_psor_residuals()
    , m_eq_residuals()
    , m_grad_residuals()
  {}

};

// SOLVER_PROFILING_DATA_H
#endif
