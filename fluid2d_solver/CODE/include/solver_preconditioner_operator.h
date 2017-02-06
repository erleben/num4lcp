#ifndef SOLVER_PRECONDITIONER_OPERATOR_H
#define SOLVER_PRECONDITIONER_OPERATOR_H

#include <solver_find_free.h>

#include <cusp/linear_operator.h>
#include <cusp/array1d.h>
#include <cusp/blas/blas.h>
#include <cusp/multiply.h>


template <typename preconditioner >
class PreconditionerOperator
  : public cusp::linear_operator< typename preconditioner::value_type
                                , typename preconditioner::memory_space
                                , typename preconditioner::index_type
                                >
{
public:

  typedef typename preconditioner::value_type    T;
  typedef typename preconditioner::memory_space  M;
  typedef          cusp::linear_operator<T,M>    parent_type;
  
  size_t                  const   m_N;
  preconditioner          const & m_P;
  cusp::array1d<bool,M>   const & m_active;
  cusp::array1d<bool,M>           m_free;
  
  PreconditionerOperator(
                             preconditioner        const & P
                           , cusp::array1d<bool,M> const & active
                           )
    : parent_type(P.num_rows, P.num_cols)
    , m_P(P)
    , m_active(active)
    , m_N(active.size())
  {
    assert(P.num_rows == P.num_cols || !"PreconditionOperator() Incompatible dimensions");
    
    m_free = find_free(m_active);
  }
  
  template <  typename V  >
  void operator()(  V const & v, V & w) const
  {
    V v_A(m_N, 0.0);
    V v_F(m_N, 0.0);

    // Evaluates Px = y, using the form
    //
    // | P_{AA}   0    | | v_A | = | w_{A} |
    // |   0    I_{FF} | | v_F |   | v_{F} |
    //
    cusp::blas::xmy(m_active, v, v_A); // x .* {A} = x_{A}
    cusp::multiply(m_P, v_A, w);       // w_{A} = v_tmp = Px_{A}
    cusp::blas::xmy(m_free, v, v_F);   // w_{F} = x_{F}
    cusp::blas::axpy(v_F,w,1);         // w = w_{F} + w_{A}
  }
};


template <typename M>
inline PreconditionerOperator<M> make_preconditioner_operator(M const & P, cusp::array1d<bool, typename M::memory_space> const & active)
{
  return PreconditionerOperator<M>(P, active);
}

// SOLVER_PRECONDITIONER_OPERATOR_H
#endif
