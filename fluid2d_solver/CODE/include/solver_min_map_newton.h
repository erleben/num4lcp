#ifndef SOLVER_MIN_MAP_NEWTON_H
#define SOLVER_MIN_MAP_NEWTON_H

#include <solver_status_codes.h>
#include <solver_project.h>
#include <solver_find_active.h>
#include <solver_compute_minimum_map.h>
#include <solver_preconditioner_operator.h>
#include <solver_schur_operator.h>
#include <solver_jacobian_operator.h>
#include <solver_jacobian_transpose_operator.h>
#include <solver_linear_system_solver.h>
#include <util_params.h>


#include <cusp/array1d.h>
#include <cusp/blas/blas.h>
#include <cusp/multiply.h>
#include <cusp/monitor.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/cr.h>
#include <cusp/krylov/gmres.h>
#include <cusp/krylov/bicg.h>
#include <cusp/krylov/bicgstab.h>

#include <limits>

namespace details
{

  template<typename V>
  inline void negate(V const & x, V & y)
  {
    thrust::transform( x.begin(), x.end(), y.begin(), thrust::negate<typename V::value_type>() );
  }
  
  
  /**
   *
   * @param f_0              Upon invokation holds the initial merit value at
   *                         x_0 (tau=0), upon return holds the value of f for
   *                         the given tau value.
   * @param Df               The directional derivative of merit function in
   *                         direction dx evaluated at x_0
   * @param dx               The serach direction (should be a descent
   *                         direction).
   *
   * @param alpha            The step reduction parameter
   * @param beta             The sufficient decrease parameter
   * @param gamma            Too small step length guard, squared acceptable
   *                          minimum step length.
   * @param max_iterations   Maximum allowed iterations to perform
   *
   * @param A                The coefficient matrix of y = A x + b
   * @param b                The rhs vector in  y = A x + b
   * @param x                The value of  x in y = A x + b
   * @param y                The lhs vector  of y = A x + b
   * @param H                The coefficient matrix of y = A x + b
   * @param tau              Upon return holds the resulting step length
   *                         value to use.
   * @param iteration        The number of iterations used to find tau
   * @param status           The exit-status from the line-search, can be
   *                         SUFFICIENT_DECRASE, MAX_LIMIT, TOO_SMALL_STEP.
   */
  template<typename M, typename V, typename VB, typename T>
  inline void projected_back_tracking_line_search(
                                                  T                    & f
                                                  , T            const & Df
                                                  , V            const & dx
                                                  , T            const & alpha
                                                  , T            const & beta
                                                  , T            const & gamma
                                                  , unsigned int const & max_iterations
                                                  , M            const & A
                                                  , V            const & b
                                                  , VB           const & is_bilateral
                                                  , V                  & x
                                                  , V                  & y
                                                  , V                  & H
                                                  , T                  & tau
                                                  , unsigned int       & iteration
                                                  , unsigned int       & status
                                                  )
  {
    assert(f > 0.0            || !"projected_back_tracking_line_search(): Objective function must be positive"    );
    assert(Df < 0.0           || !"projected_back_tracking_line_search(): directional derivative must be negative");
    assert(alpha > 0.0        || !"projected_back_tracking_line_search(): alpha must be positive"                 );
    assert(alpha < 1.0        || !"projected_back_tracking_line_search(): alpha must be less than one"            );
    assert(beta > 0.0         || !"projected_back_tracking_line_search(): beta must be positive"                  );
    assert(beta < 1.0         || !"projected_back_tracking_line_search(): beta must be less than one"             );
    assert(gamma > 0.0        || !"projected_back_tracking_line_search(): gamma must be positive"                 );
    assert(gamma < 0.001      || !"projected_back_tracking_line_search(): gamma should be smaller"                );
    assert(max_iterations > 0 || !"projected_back_tracking_line_search(): Maximum iterations must be positive"    );
    
    V tau_dx(x.size(), T(0.0));
    
    V x_0(x.size(), T(0.0));
    cusp::copy( x, x_0 );
    
    T const f_0 = f;
    
    tau = T(1.0);
    
    status = ITERATING;
    
    iteration = 0u;
    
    bool const done = false;
    while ( not done  )
    {
      // x_k = x + tau*dx
      cusp::copy( x_0, x );
      cusp::copy( dx, tau_dx );
      blas::scal( tau_dx, tau);
      blas::axpy( tau_dx, x, 1 );
      
      // x_k = max(0,x_k)   // Note: Only for unilateral constraints, for bilateral do nothing
      ::project( x, is_bilateral, x );
      
      // y_k = A x_k + b
      cusp::multiply( A, x, y );
      blas::axpy( b, y, T(1.0) );

      // H_k = min(y_k,x_k)     # Note: Only for unilateral constraints, for bilateral H_k = y_k
      ::compute_minimum_map( y, x, is_bilateral, H );
      
      // f_k = 1/2 H_k^T H_k
      f = T(0.5) * blas::dot( H, H );

      if(iteration >= max_iterations)
      {
        status = MAX_LIMIT;
        return;
      }
      
      if (f <= (f_0 + tau*beta*Df))
      {
        status = SUFFICIENT_DECREASE;
        return;
      }
      
      if (tau*tau < gamma)
      {
        status = TOO_SMALL_STEP;
        return;
      }
      
      tau *= alpha;
      
      ++iteration;
    }
  }
  
}//end namespace details


template <
    typename matrix_type
  , typename vector_type1
  , typename vector_type2
  , typename preconditioner
  >
inline __host__ void min_map_newton(  matrix_type    const & A
                                    , vector_type1   const & b
                                    , vector_type1         & x
                                    , vector_type2   const & is_bilateral
                                    , preconditioner       & P
                                    , Params         const & params
                                    , unsigned int         & newton_status
                                    , unsigned int         & newton_iteration
                                    , cusp::array1d<typename vector_type1::value_type, cusp::host_memory> & residuals
                                    )
{
  namespace blas = cusp::blas;
  
  using std::min;
  using std::max;
  using std::fabs;

  typedef typename matrix_type::value_type                 T;
  typedef          vector_type1                            V;
  typedef          vector_type2                            VB;
  typedef typename matrix_type::memory_space               M;
  typedef          JacobianOperator<matrix_type>           jacobian_operator_type;
  typedef          JacobianTransposeOperator<matrix_type>  jacobian_transpose_operator_type;
  typedef          PreconditionerOperator<preconditioner>  preconditioner_operator_type;
  typedef          ShurOperator<matrix_type, vector_type1> schur_operator_type;

  assert(A.num_cols == A.num_rows || !"min_map_newton(): incompatible dimensions");

  bool         const done = false;
  unsigned int const N    = A.num_rows;

  newton_status      = ITERATING;

  newton_iteration   = 0u;

  if( params.profiling() && params.record_convergence() )
  {
    residuals.resize(params.newton_max_iterations(), T(0.0));
  }

  V  neg_H(N);
  V  dx(N);
  V  grad_f(N);
  V  H(N);
  V  y(N);
  VB is_active( N );

  // y = Ax + b
  cusp::multiply( A, x, y );
  blas::axpy( b, y, 1 );

  // H = min(y,x)       # Only for unilateral constraints, for bilateral we use H = y
  compute_minimum_map( y, x, is_bilateral, H );

  T f_old = std::numeric_limits<T>::max();
  T f     = T(0.5) * blas::dot( H, H );

  while ( not done )
  {
    if( params.profiling() && params.record_convergence() )
    {
      residuals[newton_iteration] = f;
    }

    if (f < params.newton_absolute_tolerance() )
    {
      newton_status = ABSOLUTE_CONVERGENCE;
      if(params.verbose())
      {
        std::cout << "min_map_newton(): ABSOLUTE_CONVERGENCE" << std::endl;
      }
      return;
    }

    if ( fabs(f - f_old) < params.newton_relative_tolerance()*f_old )
    {
      newton_status = RELATIVE_CONVERGENCE;
      if(params.verbose())
      {
        std::cout << "min_map_newton(): RELATIVE CONVERGENCE" << std::endl;
      }
      return;
    }
    
    f_old = f;
    
    //--- Solve Newton Equation for search direction dx
    {
      details::negate(H, neg_H);
      
      find_active( y, x, is_bilateral, is_active );

      assert( params.newton_rho() >= 0.0 || !"min_map_newton(): rho must be non-negative");
      assert( params.newton_rho() < 1.0  || !"min_map_newton(): rho must less than one"  );

      T const necessary_accuracy = params.newton_rho() * blas::nrm2(H);
      T const desired_accuracy   = params.sub_solver_absolute_tolerance();
      
      unsigned int desired_max_iterations = params.sub_solver_max_iterations();
      unsigned int necessary_max_iterations = N / 10;

      schur_operator_type          S_op = make_schur_operator( A, neg_H, is_active );
      preconditioner_operator_type P_op = make_preconditioner_operator( P, is_active );
    
      unsigned int const sub_solver_max_iteration       = max( desired_max_iterations, necessary_max_iterations );
      T            const sub_solver_absolute_tolerance  = min( desired_accuracy, necessary_accuracy );
      T            const sub_solver_relative_tolerance  = params.sub_solver_relative_tolerance();
      unsigned int const sub_solver_choice              = params.sub_solver_choice();

      cusp::monitor<T> sub_solver_monitor(
                                          S_op.m_b
                                          , sub_solver_max_iteration
                                          , sub_solver_relative_tolerance
                                          , sub_solver_absolute_tolerance
                                          );

      linear_system_solver(S_op, S_op.m_b, dx, sub_solver_monitor, P_op, sub_solver_choice);

      //--- Test if sub solver found a suitable solution or not.
      if (not sub_solver_monitor.converged())
      {
        newton_status = SUB_SOLVER_FAILURE;
        if (params.verbose())
        {
          sub_solver_monitor.print();
        }
        return;
      }
    }
    
    jacobian_operator_type           J_op  = make_jacobian_operator( A, is_active );
    jacobian_transpose_operator_type JT_op = make_jacobian_transpose_operator( A, is_active );

    //---
    //--- First we compute the gradient of f at x.
    //---
    //---    df/dx = H^T J
    //---    grad f = (df/dx)^T = J^T H
    //---
    //--- Second we compute Df, the directional derivative of f in direction of dx
    //---
    //---  Df = (grad f)^T dx = H^T J dx  = H^T (J dx)
    //---
    JT_op(H, grad_f);
    T const Df = blas::dot( grad_f, dx );

    //--- Perform some fail-safe tests, to make sure the "world" is playing
    //--- nice with us
    {
      //--- Test whether the search direction is smaller than numerical precision
      if( blas::nrmmax( dx ) < params.newton_stagnation_tolerance() )
      {
        newton_status =  STAGNATION;
        if(params.verbose())
        {
          std::cout << "min_map_newton(): STAGNATION" << std::endl;
        }
        return;
      }
      //--- Test if the gradient is too close to zero-gradient
      if ( blas::nrm2( grad_f ) < params.newton_absolute_tolerance() )
      {
        newton_status =  LOCAL_MINIMA;
        if(params.verbose())
        {
          std::cout << "min_map_newton(): LOCAL MINIMA" << std::endl;
        }
        return;
      }
      //--- Test if the search directoin is a descent direction
      if ( Df >= 0.0 )
      {
        newton_status =  NON_DESCEND_DIRECTION;
        if(params.verbose())
        {
          std::cout << "min_map_newton(): NON_DESCEND_DIRECTION" << std::endl;
        }
        return;
      }
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
                                                   f
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
            std::cout << "min_map_newton(): Line search step length was too small" << std::endl;
          }
          newton_status = LINE_SEARCH_FAILURE;
          return;
        }
        if (line_search_status == MAX_LIMIT)
        {
          if(params.verbose())
          {
            std::cout << "min_map_newton(): Line search reached maximum iteration limit, giving up" << std::endl;
          }
          newton_status = LINE_SEARCH_FAILURE;
          return;
        }
      }
    }

    if(params.verbose())
    {
      std::cout << "min_map_newton(): iteration = " << newton_iteration << " f = " << f << std::endl;
    }
    
    ++newton_iteration;
  }
  
}

// SOLVER_MIN_MAP_NEWTON_H
#endif
