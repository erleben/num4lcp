#ifndef SOLVER_SCHUR_OPERATOR_H
#define SOLVER_SCHUR_OPERATOR_H

#include <solver_find_free.h>

#include <cusp/linear_operator.h>
#include <cusp/csr_matrix.h>
#include <cusp/array1d.h>
#include <cusp/blas/blas.h>
#include <cusp/multiply.h>



template <typename M, typename V>
class ShurOperator
  : public cusp::linear_operator<typename M::value_type, typename M::memory_space, typename M::index_type>
{

public:
  typedef typename M::value_type               T;
  typedef typename M::memory_space             MS;
  typedef typename M::index_type               I;
  typedef          cusp::linear_operator<T,MS> parent_type;
  
  size_t                  const   m_N;
  M                       const & m_A;
  V                               m_b;
  cusp::array1d<bool,MS>  const & m_active;
  cusp::array1d<bool,MS>          m_free;
  
  ShurOperator(M const & A, V const & b, cusp::array1d<bool,MS> const & active)
    : parent_type(A.num_rows, A.num_cols)
    , m_A(A)
    , m_b(b)
    , m_active(active)
    , m_N(active.size())
  {
    assert(A.num_rows == A.num_cols || !"ShurOperator() Incompatible dimensions"); // sanity check
    
    m_free = find_free(m_active);
    
    //
    // Given the partitioning of A x = b
    //
    //  | A_AA   A_AF | |x_A|    |b_A|
    //  | A_FA   A_FF | |x_F|  = |b_F|
    //
    // Then given A_FA = 0 and A_FF = I we have
    //
    //  | A_AA   A_AF | |x_A|    |b_A|
    //  |   0     I   | |x_F|  = |b_F|
    //
    // Using the Shur Complement we find
    //
    //  | A_AA    0 | |x_A|    |b_A - A_AF b_F|
    //  |   0     I | |x_F|  = |b_F           |
    //
    // Hence when solving this reduced system we must store
    // the modified right hand side
    //
    //        | b_A - A_AF b_F |
    //   b' = |  b_F           |
    //
    
    V c(m_N, 0.0);
    V d(m_N, 0.0);

    cusp::blas::xmy(m_free, m_b, c);        // b_{F} = {F} .* b; c_{F} = b_{F}; c_{A} = 0
    cusp::multiply(m_A, c, d);              // d = A*b_{F}
    cusp::blas::xmy(m_active,d, c);         // c_{A} = {A} .* d; c = {A} .* A_{AF} b_{F}; c_{F} = 0
    cusp::blas::axpy(c,m_b,-1);             // b_{A} = b_{A} - A_{AF}b_{F}
  }
  
  void operator()( V const & x, V & y) const
  {
    V x_active(m_N, 0.0);
    V b_free(m_N, 0.0);
    
    //
    // Given the partitioning of A x = b
    //
    //  | A_AA   A_AF | |x_A|    |b_A|
    //  | A_FA   A_FF | |x_F|  = |b_F|
    //
    // Then given A_FA = 0 and A_FF = I we have
    //
    //  | A_AA   A_AF | |x_A|    |b_A|
    //  |   0     I   | |x_F|  = |b_F|
    //
    // Using the Shur Complement we find
    //
    //  | A_AA    0 | |x_A|    |b_A - A_AF b_F|
    //  |   0     I | |x_F|  = |b_F           |
    //
    // Hence the result of this operator should yield
    //
    //       | A_AA x_A|
    //   y = |  b_F    |
    //
    
    cusp::blas::xmy(m_active, x, x_active);    // x_{A}
    cusp::multiply(m_A, x_active, y);          // y = A *  x_{A}
    cusp::blas::xmy(m_active, y, y);           // y_{F} = 0
    cusp::blas::xmy(m_free, m_b, b_free);      // b_{F} = {F} .* b; b_{A} = 0
    cusp::blas::axpy(b_free, y, 1);            // y_{F} = b_f
  }
  
};


template< typename M, typename V>
inline ShurOperator<M,V> make_schur_operator(M const & A, V const & b, cusp::array1d<bool, typename M::memory_space> const & active)
{
  return ShurOperator<M,V>(A,b,active);
}

// SOLVER_SCHUR_OPERATOR_H
#endif
