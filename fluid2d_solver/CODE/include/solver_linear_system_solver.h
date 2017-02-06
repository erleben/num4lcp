#ifndef SOLVER_LINEAR_SYSTEM_SOLVER_H
#define SOLVER_LINEAR_SYSTEM_SOLVER_H

#include <cusp/monitor.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/cr.h>
#include <cusp/krylov/gmres.h>
#include <cusp/krylov/bicg.h>
#include <cusp/krylov/bicgstab.h>
#include <thrust/fill.h>


template<
typename matrix_type
, typename vector_type
, typename monitor_type
, typename preconditioner_type
>
inline void linear_system_solver(
                          matrix_type & A
                          , vector_type & b
                          , vector_type & x
                          , monitor_type & monitor
                          , preconditioner_type & P
                          , unsigned int const & choice
                          )
{
  typedef typename vector_type::value_type T;
  
  thrust::fill(x.begin(), x.end(), T(0.0));

  switch (choice)
  {
    case 0:
    {
      cusp::krylov::cg(A, x, b, monitor, P);
    }
      break;
    case 1:
    {
      cusp::krylov::cr(A, x, b, monitor, P);
    }
      break;
      //      case 2:
      //      {
      //        // TODO: Currently our *BlockOperators does not support transpose.
      //        std::cerr << "invoke_sub_solver(): Unsupported feature. Unable to call subsolver BICG" << std::endl;
      //        // cusp::krylov::bicg(A, At, x, b, monitor, P, Pt);
      //      }
      //      break;
    case 3:
    {
      cusp::krylov::bicgstab(A, x, b, monitor, P);
    }
      break;
    case 4:
    {
      unsigned int restart = 10u;  // hardwired constant

      cusp::krylov::gmres(A, x, b, restart, monitor, P);
    }
      break;
    default:
      assert(false || !"linear_system_solver(): Unrecognized subsolver choice");
      break;
  }
}


// SOLVER_LINEAR_SYSTEM_SOLVER_H
#endif
