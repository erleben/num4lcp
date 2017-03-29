#ifndef SOLVER_GRADIENT_DESCENT_H
#define SOLVER_GRADIENT_DESCENT_H

#include <solver_status_codes.h>
#include "solver_min_map_newton.h"
#include <solver_find_active.h>
#include <solver_psor.h>

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

template <
          typename M
        , typename V
        , typename VB
        >
inline __host__ void gradient_descent(
                                      M              const & A
                                      , V            const & b
                                      , V                  & x
                                      , VB           const & is_bilateral
                                      , Params       const & params
                                      , unsigned int       & grad_status
                                      , unsigned int       & grad_iteration
                                      , cusp::array1d<typename V::value_type, cusp::host_memory> & residuals
                                      )
{
  typedef typename V::value_type T;
  typedef          JacobianTransposeOperator<M> jacobian_transpose_operator_type;

  using std::fabs;
  using std::max;

  unsigned int const N         = A.num_cols;
  bool         const done      = false;
  T            const too_small = 1e-6;

  assert( A.num_cols == A.num_rows                     || !"grad(): A matrix has incompatible dimentions"     );
  assert( params.grad_max_iterations() > 0u            || !"grad(): max iterations must be positive"          );
  assert( params.grad_absolute_tolerance() >= T(0.0)   || !"grad(): absolute tolerance must be non-negative"  );
  assert( params.grad_relative_tolerance() >= T(0.0)   || !"grad(): relative tolerance must be non-negative"  );
  assert( params.grad_stagnation_tolerance() >= T(0.0) || !"grad(): stagnation tolerance must be non-negative");

  grad_status    = ITERATING;
  grad_iteration = 0u;

  if (params.profiling() && params.record_convergence() )
  {
    residuals.resize(params.grad_max_iterations(), T(0.0) );
  }

  V y(N);
  V H(N);
  V A_diag(N);
  V grad_f(N);
  V dx(N);
  VB is_active(N);

  details::make_safe_diagonal_vector(A, A_diag, too_small);

  cusp::multiply(A,x,y);
  cusp::blas::axpy(b,y,1);

  compute_minimum_map( y, x, is_bilateral, H );
  T residual     = T(0.5) * cusp::blas::dot(H,H);
  T residual_old = std::numeric_limits<T>::max();

  while( not done )
  {
    if( grad_iteration >= params.grad_max_iterations() )
    {
      grad_status = MAX_LIMIT;
      if(params.verbose())
      {
        std::cout << "gradient_descent(): MAX_LIMIT test passed" << std::endl;
      }
      return;
    }

    find_active( y, x, is_bilateral, is_active );

    jacobian_transpose_operator_type JT_op = make_jacobian_transpose_operator( A, is_active );

    JT_op( H, grad_f );
    details::negate( grad_f, dx );

    T const Df = blas::dot( grad_f, dx );
    //--- Test whether the search direction is smaller than numerical precision
    if( blas::nrmmax( dx ) < params.grad_stagnation_tolerance() )
    {
      grad_status =  STAGNATION;
      if(params.verbose())
      {
        std::cout << "gradient_descent(): STAGNATION" << std::endl;
      }
      return;
    }

    //--- Test if the gradient is too close to zero-gradient
    if ( blas::nrm2( grad_f ) < params.grad_absolute_tolerance() )
    {
      grad_status =  LOCAL_MINIMA;
      if(params.verbose())
      {
        std::cout << "gradient_descent(): LOCAL MINIMA" << std::endl;
      }
      return;
    }
    //--- Now we are ready to perform a line search
    {
      T            const line_search_alpha          = params.line_search_alpha();
      T            const line_search_beta           = params.line_search_beta();
      T            const line_search_gamma          = params.line_search_gamma();
      unsigned int const line_search_max_iterations = params.line_search_max_iterations();
      unsigned int       line_search_status         = ITERATING;
      unsigned int       line_search_iteration      = 0u;
      T                  tau                        = T(1.0);

      details::projected_back_tracking_line_search(
                                                   residual
                                                   , Df
                                                   , dx
                                                   , line_search_alpha
                                                   , line_search_beta
                                                   , line_search_gamma
                                                   , line_search_max_iterations
                                                   , A
                                                   , b
                                                   , is_bilateral
                                                   , x
                                                   , y
                                                   , H
                                                   , tau
                                                   , line_search_iteration
                                                   , line_search_status
                                                   );

      //--- Make tests to make sure the line search was OK
      {
        if (line_search_status == TOO_SMALL_STEP)
        {
          if(params.verbose())
          {
            std::cout << "gradient_descent(): Line search step length was too small" << std::endl;
          }
          grad_status = LINE_SEARCH_FAILURE;
          return;
        }
        if (line_search_status == MAX_LIMIT)
        {
          if(params.verbose())
          {
            std::cout << "gradient_descent(): Line search reached maximum iteration limit, giving up" << std::endl;
          }
          grad_status = LINE_SEARCH_FAILURE;
          return;
        }
      }
    }

    cusp::multiply(A,x,y);
    cusp::blas::axpy(b,y,1);

    compute_minimum_map(y,x,is_bilateral, H);

    residual = T(0.5) * cusp::blas::dot(H,H);   // We are going to use the same residual as for the min_map_newton method

    if (params.profiling() && params.record_convergence() )
    {
      residuals[grad_iteration] = residual;
    }
    if( params.verbose() )
    {
      std::cout << "gradient_descent(): iteration "
                << grad_iteration
                << " with residual "
                << residual
                << std::endl;
    }

    if (residual <= params.grad_absolute_tolerance() )
    {
      grad_status = ABSOLUTE_CONVERGENCE;
      if(params.verbose())
      {
        std::cout << "gradient_descent(): ABSOLUTE convergence test passed" << std::endl;
      }
      return;
    }
    if( fabs(residual - residual_old) <= params.grad_relative_tolerance()*residual_old )
    {
      grad_status = RELATIVE_CONVERGENCE;
      if(params.verbose() )
      {
        std::cout << "gradient_descent(): RELATIVE convergence test passed" << std::endl;
      }
      return;
    }
    residual_old = residual;

    assert(residual_old >= T(0.0) || !"gradient_descent(): internal error");

    ++grad_iteration;
  }
}

// SOLVER_GRADIENT_DESCENT_H
#endif
