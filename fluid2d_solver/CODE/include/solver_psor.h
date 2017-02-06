#ifndef SOLVER_PSOR_H
#define SOLVER_PSOR_H

#include <solver_status_codes.h>
#include <solver_sparse_row_dot.h>
#include <util_params.h>

#include <cusp/array1d.h>
#include <cusp/multiply.h>
#include <cusp/blas/blas.h>
#include <cusp/csr_matrix.h>

#include <thrust/functional.h>

#include <cassert>
#include <cmath>
#include <algorithm>
#include <limits>


namespace details
{
  
  template<typename M, typename V>
  inline __host__ void make_safe_diagonal_vector(
                                                 M const & A
                                                 , V & d
                                                 , typename V::value_type const & too_small
                                                 )
  {
    using std::fabs;

    typedef typename V::value_type T;

    cusp::extract_diagonal(A, d);
    
    for (unsigned int i = 0u; i < d.size(); ++i)
    {
      d[i] = (fabs(d[i]) < too_small) ? T(1.0) : d[i];
    }
  }
  
  template <typename T>
  struct AbsoluteValue : public thrust::unary_function<T,T>
  {
    __host__ __device__ T operator()(T x)
    {
      using std::fabs;
      
      return fabs(x);
    }
  };

  template<typename V>
  inline __host__ void abs(V const & x, V & y)
  {
    typedef typename V::value_type T;
    
    thrust::transform(x.begin(), x.end(), y.begin(), AbsoluteValue<T>() );
  }
  
}// end namespace details


template <
          typename M
        , typename V
        , typename VB
        >
inline __host__ void psor(
                            M            const & A
                          , V            const & b
                          , V                  & x
                          , VB           const & is_bilateral
                          , Params       const & params
                          , unsigned int       & psor_status
                          , unsigned int       & psor_iteration
                          , cusp::array1d<typename V::value_type, cusp::host_memory> & residuals
                          )
{
  typedef typename V::value_type T;

  using std::fabs;
  using std::max;

  unsigned int const N         = A.num_cols;
  bool         const done      = false;
  T            const gamma     = params.psor_relaxation();
  T            const too_small = 1e-6;

  assert( A.num_cols == A.num_rows                     || !"psor(): A matrix has incompatible dimentions"     );
  assert( gamma > T(0.0)                               || !"psor(): relaxation value must be positive"        );
  assert( gamma < T(2.0)                               || !"psor(): relaxation value must less than two"      );
  assert( params.psor_max_iterations() > 0u            || !"psor(): max iterations must be positive"          );
  assert( params.psor_absolute_tolerance() >= T(0.0)   || !"psor(): absolute tolerance must be non-negative"  );
  assert( params.psor_relative_tolerance() >= T(0.0)   || !"psor(): relative tolerance must be non-negative"  );
  assert( params.psor_stagnation_tolerance() >= T(0.0) || !"psor(): stagnation tolerance must be non-negative");

  psor_status    = ITERATING;
  psor_iteration = 0u;

  if (params.profiling() && params.record_convergence() )
  {
    residuals.resize(params.psor_max_iterations(), T(0.0) );
  }

  V y(N);
  V H(N);
  V A_diag(N);

  details::make_safe_diagonal_vector(A, A_diag, too_small);

  T residual_old = std::numeric_limits<T>::max();

  T dx = T(0.0);    // Largest absolute change in any x-coordinate, used to test for stagnation
  
  while( not done )
  {
    if( psor_iteration >= params.psor_max_iterations() )
    {
      psor_status = MAX_LIMIT;
      if(params.verbose())
      {
        std::cout << "psor(): MAX_LIMIT test passed" << std::endl;
      }
      return;
    }

    dx = T(0.0);

    for (unsigned int i = 0u; i < N; ++i)
    {
      T const old_xi = x[i];
      
      T row_product = T(0.0);
      sparse_row_dot(A, x, i, row_product);
      
      T const ri = b[i] + row_product;
      
      T const update_value = x[i] -  gamma*ri / A_diag[i] ;

      x[i] = is_bilateral[i] ? update_value : max( T(0.0), update_value );
      
      dx = max(dx, fabs(x[i] - old_xi));
    }

    cusp::multiply(A,x,y);
    cusp::blas::axpy(b,y,1);

    compute_minimum_map(y,x,is_bilateral, H);

    T const residual = T(0.5) * cusp::blas::dot(H,H);   // We are going to use the same residual as for the min_map_newton method

    if (params.profiling() && params.record_convergence() )
    {
      residuals[psor_iteration] = residual;
    }
    if( params.verbose() )
    {
      std::cout << "psor(): iteration "
                << psor_iteration
                << " with residual "
                << residual
                << std::endl;
    }

    if (residual <= params.psor_absolute_tolerance() )
    {
      psor_status = ABSOLUTE_CONVERGENCE;
      if(params.verbose())
      {
        std::cout << "psor(): ABSOLUTE convergence test passed" << std::endl;
      }
      return;
    }
    if( fabs(residual - residual_old) <= params.psor_relative_tolerance()*residual_old )
    {
      psor_status = RELATIVE_CONVERGENCE;
      if(params.verbose() )
      {
        std::cout << "psor(): RELATIVE convergence test passed" << std::endl;
      }
      return;
    }
    if (dx < params.psor_stagnation_tolerance() )
    {
      psor_status = STAGNATION;
      if(params.verbose())
      {
        std::cout << "psor(): STAGNATION convergence test passed" << std::endl;
      }
      break;
    }
    residual_old = residual;

    assert(residual_old >= T(0.0) || !"psor(): internal error");

    ++psor_iteration;
  }
}

// SOLVER_PSOR_H
#endif
