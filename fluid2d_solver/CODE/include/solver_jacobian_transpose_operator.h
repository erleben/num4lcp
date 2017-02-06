#ifndef SOLVER_JACOBIAN_TRANSPOSE_OPERATOR_H
#define SOLVER_JACOBIAN_TRANSPOSE_OPERATOR_H

#include <solver_find_free.h>

#include <cusp/linear_operator.h>
#include <cusp/array1d.h>
#include <cusp/blas/blas.h>
#include <cusp/multiply.h>


/**
 * This functor assumes that the A-matrix is symmetric
 *
 */
template <typename Matrix>
class JacobianTransposeOperator
  : public cusp::linear_operator<
                                  typename Matrix::value_type
                                , typename Matrix::memory_space
                                , typename Matrix::index_type
                                >
{
public:
  
  typedef typename Matrix::value_type   T;
  typedef typename Matrix::memory_space M;
  typedef typename Matrix::index_type   I;
  typedef          cusp::linear_operator<T, M, I> parent_type;
  
  Matrix                   const & m_A;
  cusp::array1d<bool, M>   const & m_active;
  cusp::array1d<bool, M>           m_free;
  
public:
  
  JacobianTransposeOperator(
                     Matrix                const & A
                   , cusp::array1d<bool,M> const & active_mask
                   )
    : parent_type(A.num_rows, A.num_cols)
    , m_A(A)
    , m_active(active_mask)
  {
    assert(A.num_rows == A.num_cols || !"JacobianTransposeOperator(): Only square matrices are supported");
    
    m_free = find_free(m_active);
  }
  
  template <typename V >
  void operator()(  V const & x, V & y) const
  {
    // We will compute
    //
    // y = J^T x
    //
    //  |y_A|      |A_AA  A_AF |^T | x_A |
    //  |y_F|  =   | 0    I    |   | x_F |
    //
    //  |y_A|      | A_AA   0 | | x_A |
    //  |y_F|  =   | A_FA   I | | x_F |
    //
    //   y_A  =     A_AA x_A
    //
    //                           | x_A |
    //  |y_F|  =   | A_FA  I |   | x_F |  =   A_FA x_A + x_F
    //
    //  |y_F|  =   A_FA x_A + x_F
    //
    //
    V tmp(x.size());
    cusp::blas::xmy(m_active, x, y);
    cusp::multiply(m_A, y, tmp);       // tmp = A * [x_A, 0]
    cusp::blas::xmy(m_free, x, y);     // y   = [ 0 , x_free ]
    cusp::blas::axpy(tmp, y, 1);       // y   = [ A_AA * x_active,   A_FA * x_active +  x_free ]
  }
};


template <typename M>
inline JacobianTransposeOperator<M> make_jacobian_transpose_operator(M const & A, cusp::array1d<bool,typename M::memory_space> const & active)
{
  return JacobianTransposeOperator<M>(A, active);
}


// SOLVER_JACOBIAN_TRANSPOSE_OPERATOR_H
#endif
