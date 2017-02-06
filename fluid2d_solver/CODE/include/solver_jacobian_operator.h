#ifndef SOLVER_JACOBIAN_OPERATOR_H
#define SOLVER_JACOBIAN_OPERATOR_H

#include <solver_find_free.h>

#include <cusp/linear_operator.h>
#include <cusp/array1d.h>
#include <cusp/blas/blas.h>
#include <cusp/multiply.h>


template <typename Matrix>
class JacobianOperator
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
  
  JacobianOperator(
                     Matrix                const & A
                   , cusp::array1d<bool,M> const & active_mask
                   )
    : parent_type(A.num_rows, A.num_cols)
    , m_A(A)
    , m_active(active_mask)
  {
    assert(A.num_rows == A.num_cols || !"JacobianOperator(): Only square matrices are supported");
    
    m_free = find_free(m_active);
  }
  
  template <typename V >
  void operator()(  V const & dx
                  , V       & Jdx
                  ) const
  {
    V v_tmp(dx.size(), 0.0); 
                                      
    cusp::multiply(m_A, dx, Jdx);        // Adx = A*dx
    cusp::blas::xmy(m_active, Jdx, Jdx); // Jdx = active .* Jdx | Jdx_{F} = 0
    cusp::blas::xmy(m_free, dx, v_tmp);  // v_tmp = dx_F <- free .* dx
    cusp::blas::axpy(v_tmp, Jdx, 1);     // Jdx = Jdx + dx_{F}
  }
};


template <typename M>
inline JacobianOperator<M> make_jacobian_operator(M const & A, cusp::array1d<bool,typename M::memory_space> const & active)
{
  return JacobianOperator<M>(A, active);
}


// SOLVER_JACOBIAN_OPERATOR_H
#endif
