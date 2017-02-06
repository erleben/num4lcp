#ifndef UNIT_TESTING_H
#define UNIT_TESTING_H

#include <util_params.h>

#include <solver_compute_minimum_map.h>
#include <solver_find_active.h>
#include <solver_find_free.h>
#include <solver_project.h>
#include <solver_jacobian_operator.h>
#include <solver_jacobian_transpose_operator.h>
#include <solver_preconditioner_operator.h>
#include <solver_sparse_row_dot.h>
#include <solver_psor.h>
#include <solver_schur_operator.h>
#include <solver_min_map_newton.h>
#include <solver_solve_lcp.h>
#include <solver_solve_eq.h>

#include <cusp/gallery/poisson.h>
#include <cusp/print.h>
#include <cusp/linear_operator.h>
#include <cusp/csr_matrix.h>
#include <cusp/precond/ainv.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>


template<typename T>
inline bool CHECK_CLOSE(T const & a, T const & b, T const & eps, std::string const & message = "")
{
  using std::fabs;

  T const too_small = 1e-5;

  if(fabs(b) < too_small &&  fabs(a)< too_small)
    return true;
  else if( fabs(a - b)<=eps*fabs(a) && fabs(a - b)<= eps*fabs(b) )
    return true;

  std::cout << "CHECK_CLOSE(" << a << "," << b << ") failed : " << message << std::endl;

  return false;
}




inline void test_minimum_map()
{
  typedef float T;

  std::cout << "############# testing minimum map #################" << std::endl;
  
  typedef cusp::host_memory   M;
  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> y(5);
  cusp::array1d<float, M> H(5);
  
  x[0] = -1.0;
  x[1] =  1.0;
  x[2] = -1.0;
  x[3] =  1.0;
  x[4] = -1.0;

  y[0] =  1.0;
  y[1] = -1.0;
  y[2] =  1.0;
  y[3] = -1.0;
  y[4] =  1.0;

  H[0] =  0.0;
  H[1] =  0.0;
  H[2] =  0.0;
  H[3] =  0.0;
  H[4] =  0.0;

  cusp::array1d<bool, M> is_bilateral(5, false);

  compute_minimum_map(y, x, is_bilateral, H);
  
  std::cout << "Result of H = min(x,y) = (-1,-1,-1,-1,-1):" << std::endl;
  CHECK_CLOSE( H[0], T(-1.0), T(0.1), "Minimum map test 1: H[0] failed");
  CHECK_CLOSE( H[1], T(-1.0), T(0.1), "Minimum map test 1: H[1] failed");
  CHECK_CLOSE( H[2], T(-1.0), T(0.1), "Minimum map test 1: H[2] failed");
  CHECK_CLOSE( H[3], T(-1.0), T(0.1), "Minimum map test 1: H[3] failed");
  CHECK_CLOSE( H[4], T(-1.0), T(0.1), "Minimum map test 1: H[4] failed");

  is_bilateral[0] = true;
  is_bilateral[1] = true;
  is_bilateral[2] = true;
  is_bilateral[3] = false;
  is_bilateral[4] = false;

  compute_minimum_map(y, x, is_bilateral, H);

  std::cout << "Result of H = min(x,y) = ( 1,-1, 1,-1,-1):" << std::endl;
  CHECK_CLOSE( H[0], T( 1.0), T(0.1), "Minimum map test 1: H[0] failed");
  CHECK_CLOSE( H[1], T(-1.0), T(0.1), "Minimum map test 1: H[1] failed");
  CHECK_CLOSE( H[2], T( 1.0), T(0.1), "Minimum map test 1: H[2] failed");
  CHECK_CLOSE( H[3], T(-1.0), T(0.1), "Minimum map test 1: H[3] failed");
  CHECK_CLOSE( H[4], T(-1.0), T(0.1), "Minimum map test 1: H[4] failed");
}

inline void test_find_active()
{
  std::cout << "############# testing find active #################" << std::endl;
  
  typedef cusp::host_memory M;
  
  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> y(5);
  cusp::array1d<bool, M> is_solid(5);
  cusp::array1d<bool, M> is_active(5);
  cusp::array1d<bool, M> is_free(5);
  
  x[0] = -1.0;   y[0] =  1.0;  is_solid[0] =  true;
  x[1] =  1.0;   y[1] = -1.0;  is_solid[1] =  true;
  x[2] = -1.0;   y[2] =  1.0;  is_solid[2] =  false;
  x[3] =  1.0;   y[3] = -1.0;  is_solid[3] =  false;
  x[4] = -1.0;   y[4] =  1.0;  is_solid[4] =  false;
  
  find_active(y, x, is_solid, is_active);
  
  is_free = find_free(is_active);

  std::cout << "is_active: (true, true, false, true, false)" << std::endl;
  cusp::print(is_active);
  
  std::cout << "is_free: (false, false, true, false, true)" << std::endl;
  cusp::print(is_free);
}

inline void test_project()
{
  std::cout << "############# testing project #################" << std::endl;
  
  typedef cusp::host_memory M;
  
  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> y(5);
  cusp::array1d<bool, M> is_bilateral(5, false);

  x[0] = -1.0;
  x[1] =  1.0;
  x[2] = -1.0;
  x[3] =  1.0;
  x[4] = -1.0;
  
  project(x,is_bilateral, y);

  std::cout << "Result of projection, y = max(0,x) should be (0,1,0,1,0):" << std::endl;
  cusp::print(y);

  x[0] = -1.0;  is_bilateral[0] = true;
  x[1] =  1.0;  is_bilateral[1] = true;
  x[2] = -1.0;  is_bilateral[2] = true;
  x[3] =  1.0;  is_bilateral[3] = false;
  x[4] = -1.0;  is_bilateral[4] = false;
  project(x,is_bilateral, y);

  std::cout << "Result of projection, y = max(0,x) should be (-1,1,-1,1,0):" << std::endl;
  cusp::print(y);
}

inline void test_jacobian_operator()
{
  std::cout << "############# testing jacobian operator #################" << std::endl;
  
  typedef cusp::host_memory M;
  typedef JacobianOperator< cusp::csr_matrix<int, float, M> > jacobian_type;
  
  cusp::csr_matrix<int, float, M> A(5,5,13);
  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> y(5);
  cusp::array1d<bool, M> is_active(5);
  
  A.row_offsets[0] = 0;  // first offset is always zero
  A.row_offsets[1] = 2;
  A.row_offsets[2] = 5;
  A.row_offsets[3] = 8;
  A.row_offsets[4] = 11;
  A.row_offsets[5] = 13; // last offset is always num_entries
  
  A.column_indices[0] = 0;  A.values[0] = 2;
  A.column_indices[1] = 1;  A.values[1] = 1;
  A.column_indices[2] = 0;  A.values[2] = 1;
  A.column_indices[3] = 1;  A.values[3] = 2;
  A.column_indices[4] = 2;  A.values[4] = 1;
  A.column_indices[5] = 1;  A.values[5] = 1;
  A.column_indices[6] = 2;  A.values[6] = 2;
  A.column_indices[7] = 3;  A.values[7] = 1;
  A.column_indices[8] = 2;  A.values[8] = 1;
  A.column_indices[9] = 3;  A.values[9] = 2;
  A.column_indices[10] = 4; A.values[10] = 1;
  A.column_indices[11] = 3; A.values[11] = 1;
  A.column_indices[12] = 4; A.values[12] = 2;
  
  x[0] =  1.0;
  x[1] =  1.0;
  x[2] =  1.0;
  x[3] =  1.0;
  x[4] =  1.0;
  is_active[0] =  true;
  is_active[1] =  true;
  is_active[2] =  true;
  is_active[3] =  false;
  is_active[4] =  false;
  
  jacobian_type J = make_jacobian_operator(A, is_active);
  J(x, y);

  std::cout << "Result of Jacobian Operator, y = J*x, should be (3,4,4,1,1) " << std::endl;
  cusp::print(y);
}

inline void test_jacobian_transpose_operator()
{
  std::cout << "############# testing jacobian transpose operator #################" << std::endl;

  typedef cusp::host_memory M;
  typedef JacobianTransposeOperator< cusp::csr_matrix<int, float, M> > jacobian_transpose_type;

  cusp::csr_matrix<int, float, M> A(5,5,13);
  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> y(5);
  cusp::array1d<bool, M> is_active(5);

  A.row_offsets[0] = 0;  // first offset is always zero
  A.row_offsets[1] = 2;
  A.row_offsets[2] = 5;
  A.row_offsets[3] = 8;
  A.row_offsets[4] = 11;
  A.row_offsets[5] = 13; // last offset is always num_entries

  A.column_indices[0] = 0;  A.values[0] = 2;
  A.column_indices[1] = 1;  A.values[1] = 1;
  A.column_indices[2] = 0;  A.values[2] = 1;
  A.column_indices[3] = 1;  A.values[3] = 2;
  A.column_indices[4] = 2;  A.values[4] = 1;
  A.column_indices[5] = 1;  A.values[5] = 1;
  A.column_indices[6] = 2;  A.values[6] = 2;
  A.column_indices[7] = 3;  A.values[7] = 1;
  A.column_indices[8] = 2;  A.values[8] = 1;
  A.column_indices[9] = 3;  A.values[9] = 2;
  A.column_indices[10] = 4; A.values[10] = 1;
  A.column_indices[11] = 3; A.values[11] = 1;
  A.column_indices[12] = 4; A.values[12] = 2;

  x[0] =  1.0;
  x[1] =  1.0;
  x[2] =  1.0;
  x[3] =  1.0;
  x[4] =  1.0;
  is_active[0] =  true;
  is_active[1] =  true;
  is_active[2] =  true;
  is_active[3] =  false;
  is_active[4] =  false;

  jacobian_transpose_type JT = make_jacobian_transpose_operator(A, is_active);
  JT(x, y);

//          2 1 0 0 0
// A =      1 2 1 0 0
//          0 1 2 1 0
//          0 0 1 2 1
//          0 0 0 1 2
//
//          2 1 0 0 0
// J =      1 2 1 0 0
//          0 1 2 1 0
//          0 0 0 1 0
//          0 0 0 0 1
//
//          2 1 0 0 0
// JT =     1 2 1 0 0
//          0 1 2 0 0
//          0 0 1 1 0
//          0 0 0 0 1
//
//          1            3
//          1            4
//   JT *   1     =      3
//          1            2
//          1            1

  std::cout << "Result of Jacobian Transpose operator, y = JT * x should be (3,4,3,2,1) " << std::endl;
  cusp::print(y);
}

inline void test_preconditioner_operator()
{
  std::cout << "############# testing preconditioner operator #################" << std::endl;
  
  typedef cusp::host_memory M;
  typedef cusp::precond::diagonal<float, M> preconditioner_type;
  typedef PreconditionerOperator< preconditioner_type > preconditioner_operator_type;
  
  cusp::csr_matrix<int, float, M> A(5,5,13);
  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> y(5);
  cusp::array1d<bool, M> is_active(5);
  
  A.row_offsets[0] = 0;  // first offset is always zero
  A.row_offsets[1] = 2;
  A.row_offsets[2] = 5;
  A.row_offsets[3] = 8;
  A.row_offsets[4] = 11;
  A.row_offsets[5] = 13; // last offset is always num_entries
  
  A.column_indices[0] = 0;  A.values[0] = 2;
  A.column_indices[1] = 1;  A.values[1] = 1;
  A.column_indices[2] = 0;  A.values[2] = 1;
  A.column_indices[3] = 1;  A.values[3] = 2;
  A.column_indices[4] = 2;  A.values[4] = 1;
  A.column_indices[5] = 1;  A.values[5] = 1;
  A.column_indices[6] = 2;  A.values[6] = 2;
  A.column_indices[7] = 3;  A.values[7] = 1;
  A.column_indices[8] = 2;  A.values[8] = 1;
  A.column_indices[9] = 3;  A.values[9] = 2;
  A.column_indices[10] = 4; A.values[10] = 1;
  A.column_indices[11] = 3; A.values[11] = 1;
  A.column_indices[12] = 4; A.values[12] = 2;
  
  x[0] =  1.0;
  x[1] =  1.0;
  x[2] =  1.0;
  x[3] =  1.0;
  x[4] =  1.0;
  is_active[0] =  true;
  is_active[1] =  true;
  is_active[2] =  true;
  is_active[3] =  false;
  is_active[4] =  false;
  
  preconditioner_type P(A);
  preconditioner_operator_type P_op = make_preconditioner_operator(P, is_active);
  P_op(x, y);

  std::cout << "Result of applying preconditioner should be ,y = P*x,  (0.5, 0.5, 0.5, 1, 1)" << std::endl;
  cusp::print(y);
}

inline void test_sparse_row()
{
  std::cout << "############# testing sparse row dot #################" << std::endl;
  
  typedef cusp::host_memory M;
  
  cusp::csr_matrix<int, float, M> A(5,5,13);
  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> y(5);
  
  A.row_offsets[0] = 0;  // first offset is always zero
  A.row_offsets[1] = 2;
  A.row_offsets[2] = 5;
  A.row_offsets[3] = 8;
  A.row_offsets[4] = 11;
  A.row_offsets[5] = 13; // last offset is always num_entries
  
  A.column_indices[0] = 0;  A.values[0] = 2;
  A.column_indices[1] = 1;  A.values[1] = 1;
  A.column_indices[2] = 0;  A.values[2] = 1;
  A.column_indices[3] = 1;  A.values[3] = 2;
  A.column_indices[4] = 2;  A.values[4] = 1;
  A.column_indices[5] = 1;  A.values[5] = 1;
  A.column_indices[6] = 2;  A.values[6] = 2;
  A.column_indices[7] = 3;  A.values[7] = 1;
  A.column_indices[8] = 2;  A.values[8] = 1;
  A.column_indices[9] = 3;  A.values[9] = 2;
  A.column_indices[10] = 4; A.values[10] = 1;
  A.column_indices[11] = 3; A.values[11] = 1;
  A.column_indices[12] = 4; A.values[12] = 2;
  
  x[0] =  1.0;
  x[1] =  1.0;
  x[2] =  1.0;
  x[3] =  1.0;
  x[4] =  1.0;
  
  std::cout << "Solution of all sparse row products should yield (3,4,4,4,3)" << std::endl;
  for(int i=0; i<5; ++i)
  {
    float value = 0.0;
    sparse_row_dot(A,x,i, value);
    std::cout << "sum_j A_ij^T * x_j = " << value << std::endl;
  }
}

inline void test_psor()
{
  Params params;

  params.verbose() = false;
  params.psor_max_iterations() = 100;

  std::cout << "############# testing psor #################" << std::endl;
  
  typedef cusp::host_memory M;
  
  cusp::csr_matrix<int, float, M> A(5,5,13);

  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> y(5);
  cusp::array1d<float, M> b(5);

  cusp::array1d<bool, M> is_bilateral(5);
  
  is_bilateral[0] = false;
  is_bilateral[1] = false;
  is_bilateral[2] = false;
  is_bilateral[3] = false;
  is_bilateral[4] = false;
  
  A.row_offsets[0] = 0;  // first offset is always zero
  A.row_offsets[1] = 2;
  A.row_offsets[2] = 5;
  A.row_offsets[3] = 8;
  A.row_offsets[4] = 11;
  A.row_offsets[5] = 13; // last offset is always num_entries
  
  A.column_indices[0] = 0;  A.values[0] = 2;
  A.column_indices[1] = 1;  A.values[1] = 1;
  A.column_indices[2] = 0;  A.values[2] = 1;
  A.column_indices[3] = 1;  A.values[3] = 2;
  A.column_indices[4] = 2;  A.values[4] = 1;
  A.column_indices[5] = 1;  A.values[5] = 1;
  A.column_indices[6] = 2;  A.values[6] = 2;
  A.column_indices[7] = 3;  A.values[7] = 1;
  A.column_indices[8] = 2;  A.values[8] = 1;
  A.column_indices[9] = 3;  A.values[9] = 2;
  A.column_indices[10] = 4; A.values[10] = 1;
  A.column_indices[11] = 3; A.values[11] = 1;
  A.column_indices[12] = 4; A.values[12] = 2;
  
  x[0] =  1.0;
  x[1] =  0.0;
  x[2] =  1.0;
  x[3] =  0.0;
  x[4] =  1.0;

  b[0] =  0.0;
  b[1] =  1.0;
  b[2] =  0.0;
  b[3] =  1.0;
  b[4] =  0.0;
  
  cusp::multiply(A, x, y);
  cusp::blas::axpy(y,b,-1);

  x[0] =  0.0;
  x[1] =  0.0;
  x[2] =  0.0;
  x[3] =  0.0;
  x[4] =  0.0;
  {
    unsigned int psor_status    = 0u;
    unsigned int psor_iteration = 0u;
    cusp::array1d<float, cusp::host_memory> residuals;

    psor(A, b, x, is_bilateral, params, psor_status, psor_iteration, residuals);

    std::cout << "LCP solution x should be (1,0,1,0,1)" << std::endl;
    cusp::print(x);
  }
  
  is_bilateral[0] = true;
  is_bilateral[1] = true;
  is_bilateral[2] = true;
  is_bilateral[3] = true;
  is_bilateral[4] = true;
  
  {
    unsigned int psor_status    = 0u;
    unsigned int psor_iteration = 0u;
    cusp::array1d<float, cusp::host_memory> residuals;

    psor(A, b, x, is_bilateral, params, psor_status, psor_iteration, residuals);

    std::cout << "Solution x for A x + b = 0 should be (2,-2,3,-2,2)" << std::endl;
    cusp::print(x);
  }
  
}

inline void test_schur_operator()
{
  std::cout << "############# testing schur operator #################" << std::endl;
  
  typedef cusp::host_memory M;
  
  typedef ShurOperator< cusp::csr_matrix<int, float, M>, cusp::array1d<float, M> > schur_operator_type;
  
  cusp::csr_matrix<int, float, M> A(5,5,13);
  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> b(5);
  cusp::array1d<float, M> y(5);
  cusp::array1d<bool, M> is_active(5);
  
  A.row_offsets[0] = 0;  // first offset is always zero
  A.row_offsets[1] = 2;
  A.row_offsets[2] = 5;
  A.row_offsets[3] = 8;
  A.row_offsets[4] = 11;
  A.row_offsets[5] = 13; // last offset is always num_entries
  
  A.column_indices[0] = 0;  A.values[0] = 2;
  A.column_indices[1] = 1;  A.values[1] = 1;
  A.column_indices[2] = 0;  A.values[2] = 1;
  A.column_indices[3] = 1;  A.values[3] = 2;
  A.column_indices[4] = 2;  A.values[4] = 1;
  A.column_indices[5] = 1;  A.values[5] = 1;
  A.column_indices[6] = 2;  A.values[6] = 2;
  A.column_indices[7] = 3;  A.values[7] = 1;
  A.column_indices[8] = 2;  A.values[8] = 1;
  A.column_indices[9] = 3;  A.values[9] = 2;
  A.column_indices[10] = 4; A.values[10] = 1;
  A.column_indices[11] = 3; A.values[11] = 1;
  A.column_indices[12] = 4; A.values[12] = 2;
  
  b[0] =  1.0;
  b[1] =  1.0;
  b[2] =  1.0;
  b[3] =  1.0;
  b[4] =  1.0;
  
  is_active[0] =  true;
  is_active[1] =  true;
  is_active[2] =  true;
  is_active[3] =  false;
  is_active[4] =  false;
  
  schur_operator_type S_op = make_schur_operator(A, b, is_active);

  std::cout << "New modified b of Shur operator should be (1,1,0,1,1):" << std::endl;
  cusp::print(S_op.m_b);

  x[0] =  1.0;
  x[1] =  1.0;
  x[2] =  1.0;
  x[3] =  1.0;
  x[4] =  1.0;

  S_op(x, y);
  
  std::cout << "Shur operator, S*x, should yield  (3,4,3,1,1):" << std::endl;
  cusp::print(y);
}

inline void test_params()
{
  Params params;
  
  std::cout << print_params(params);
  
  params.load("default.cfg");
  
  
  std::cout << print_params(params);
}

inline void test_min_map()
{

  std::cout << "############# testing min map #################" << std::endl;

  Params params;

  typedef cusp::host_memory M;
  typedef cusp::identity_operator<float, M, int > preconditioner_type;
  
  cusp::csr_matrix<int, float, M> A(5,5,13);
  cusp::array1d<float, M> x(5);
  cusp::array1d<float, M> y(5);
  cusp::array1d<float, M> b(5);
  cusp::array1d<bool, M> is_bilateral(5);
  
  is_bilateral[0] = false;
  is_bilateral[1] = false;
  is_bilateral[2] = false;
  is_bilateral[3] = false;
  is_bilateral[4] = false;
  
  A.row_offsets[0] = 0;  // first offset is always zero
  A.row_offsets[1] = 2;
  A.row_offsets[2] = 5;
  A.row_offsets[3] = 8;
  A.row_offsets[4] = 11;
  A.row_offsets[5] = 13; // last offset is always num_entries
  
  A.column_indices[0] = 0;  A.values[0] = 2;
  A.column_indices[1] = 1;  A.values[1] = 1;
  A.column_indices[2] = 0;  A.values[2] = 1;
  A.column_indices[3] = 1;  A.values[3] = 2;
  A.column_indices[4] = 2;  A.values[4] = 1;
  A.column_indices[5] = 1;  A.values[5] = 1;
  A.column_indices[6] = 2;  A.values[6] = 2;
  A.column_indices[7] = 3;  A.values[7] = 1;
  A.column_indices[8] = 2;  A.values[8] = 1;
  A.column_indices[9] = 3;  A.values[9] = 2;
  A.column_indices[10] = 4; A.values[10] = 1;
  A.column_indices[11] = 3; A.values[11] = 1;
  A.column_indices[12] = 4; A.values[12] = 2;
  
  x[0] =  1.0;
  x[1] =  0.0;
  x[2] =  1.0;
  x[3] =  0.0;
  x[4] =  1.0;
  
  b[0] =  0.0;
  b[1] =  1.0;
  b[2] =  0.0;
  b[3] =  1.0;
  b[4] =  0.0;
  
  cusp::multiply(A, x, y);
  cusp::blas::axpy(y,b,-1);

  x[0] =  0.0;
  x[1] =  0.0;
  x[2] =  0.0;
  x[3] =  0.0;
  x[4] =  0.0;
  {
    preconditioner_type P(5,5);

    unsigned int newton_status    = 0u;
    unsigned int newton_iteration = 0u;
    cusp::array1d<float, M> residuals;

    min_map_newton(A, b, x, is_bilateral, P, params, newton_status, newton_iteration, residuals);

    std::cout << "LCP solution x should be (1,0,1,0,1)" << std::endl;
    cusp::print(x);
  }

  x[0] =  0.0;
  x[1] =  0.0;
  x[2] =  0.0;
  x[3] =  0.0;
  x[4] =  0.0;
  is_bilateral[0] = true;
  is_bilateral[1] = true;
  is_bilateral[2] = true;
  is_bilateral[3] = true;
  is_bilateral[4] = true;
  
  {
    preconditioner_type P(5,5);

    unsigned int newton_status    = 0u;
    unsigned int newton_iteration = 0u;
    cusp::array1d<float, M> residuals;

    min_map_newton(A, b, x, is_bilateral, P, params, newton_status, newton_iteration, residuals);
    
    std::cout << "Solution x for A x + b = 0 should be (2,-2, 3,-2,2)" << std::endl;
    cusp::print(x);
  }
  
}

inline void test_larger_min_map()
{
  typedef float             T;
  typedef cusp::host_memory M;
  typedef int               I;

  Params params;

  params.verbose();
  params.newton_max_iterations() = 100;

  typedef cusp::identity_operator<T, M, I > preconditioner_type;

  unsigned int grid_size[] = {32, 64, 128, 256};

  for (unsigned int p=0;p<4u;++p)
  {
    cusp::csr_matrix<I, T, M> A;

    cusp::gallery::poisson5pt(A, grid_size[p], grid_size[p]);

    unsigned int const N = A.num_rows;

    params.sub_solver_max_iterations() = N;

    cusp::array1d<T, M> x(N, T(0.0));
    cusp::array1d<T, M> y(N, T(0.0));
    cusp::array1d<T, M> b(N, T(0.0));
    cusp::array1d<bool, M> is_bilateral(N, false);

    for (unsigned int i=0u; i < N;++i)
      if(i%2==0)
        x[i] = T(1.0*i/N);
      else
        y[i] = T(1.0*i/N);

    cusp::multiply(A, x, y);
    cusp::blas::axpy(y,b,-1);

    preconditioner_type P(N,N);

    unsigned int newton_status    = 0u;
    unsigned int newton_iteration = 0u;
    cusp::array1d<T, M> residuals;

    cusp::array1d<T, M> z(N, T(0.0));
    min_map_newton(A, b, z, is_bilateral, P, params, newton_status, newton_iteration, residuals);

    for (unsigned int i=0u; i < N;++i)
    {
      CHECK_CLOSE(z[i], x[i], T(0.01), "Large test failed");
    }

  }
}

inline void test_solve_lcp()
{

  std::cout << "############# testing solve_lcp #################" << std::endl;

  Params params;

  std::vector<unsigned int> row(13u);
  std::vector<unsigned int> column(13u);
  std::vector<float> values(13u);

  row[0] = 0;  column[0] = 0;  values[0] = 2;
  row[1] = 0;  column[1] = 1;  values[1] = 1;
  row[2] = 1;  column[2] = 0;  values[2] = 1;
  row[3] = 1;  column[3] = 1;  values[3] = 2;
  row[4] = 1;  column[4] = 2;  values[4] = 1;
  row[5] = 2;  column[5] = 1;  values[5] = 1;
  row[6] = 2;  column[6] = 2;  values[6] = 2;
  row[7] = 2;  column[7] = 3;  values[7] = 1;
  row[8] = 3;  column[8] = 2;  values[8] = 1;
  row[9] = 3;  column[9] = 3;  values[9] = 2;
  row[10] = 3;  column[10] = 4; values[10] = 1;
  row[11] = 4;  column[11] = 3; values[11] = 1;
  row[12] = 4;  column[12] = 4; values[12] = 2;

  std::vector<bool> is_bilateral(5);
  is_bilateral[0] = false;
  is_bilateral[1] = false;
  is_bilateral[2] = false;
  is_bilateral[3] = false;
  is_bilateral[4] = false;

  std::vector<float> b(5);
  b[0] =   2.0;
  b[1] =   1.0;
  b[2] =   2.0;
  b[3] =   1.0;
  b[4] =   2.0;

  std::vector<float> x(5);
  x[0] =  0.0;
  x[1] =  0.0;
  x[2] =  0.0;
  x[3] =  0.0;
  x[4] =  0.0;

  ProfilingData profiling_data;

  solve_lcp_host(
                 5
                 , 13
                 , row
                 , column
                 , values
                 , x
                 , b
                 , is_bilateral
                 , params
                 , profiling_data
                 );

  std::cout << "x = [" << x[0];

  for(unsigned int i = 1u; i< 5u;++i)
    std::cout << " " << x[i];

  std::cout << " ] " << std::endl;

  x[0] =  0.0;
  x[1] =  0.0;
  x[2] =  0.0;
  x[3] =  0.0;
  x[4] =  0.0;

  solve_lcp_device(
                 5
                 , 13
                 , row
                 , column
                 , values
                 , x
                 , b
                 , is_bilateral
                 , params
                 , profiling_data
                 );

  std::cout << "x = [" << x[0];

  for(unsigned int i = 1u; i< 5u;++i)
    std::cout << " " << x[i];

  std::cout << " ] " << std::endl;
}


inline void test_solve_eq()
{

  std::cout << "############# testing solve_eq #################" << std::endl;

  Params params;

  std::vector<unsigned int> row(13u);
  std::vector<unsigned int> column(13u);
  std::vector<float> values(13u);

  row[0] = 0;  column[0] = 0;  values[0] = 2;
  row[1] = 0;  column[1] = 1;  values[1] = 1;
  row[2] = 1;  column[2] = 0;  values[2] = 1;
  row[3] = 1;  column[3] = 1;  values[3] = 2;
  row[4] = 1;  column[4] = 2;  values[4] = 1;
  row[5] = 2;  column[5] = 1;  values[5] = 1;
  row[6] = 2;  column[6] = 2;  values[6] = 2;
  row[7] = 2;  column[7] = 3;  values[7] = 1;
  row[8] = 3;  column[8] = 2;  values[8] = 1;
  row[9] = 3;  column[9] = 3;  values[9] = 2;
  row[10] = 3;  column[10] = 4; values[10] = 1;
  row[11] = 4;  column[11] = 3; values[11] = 1;
  row[12] = 4;  column[12] = 4; values[12] = 2;

  std::vector<bool> is_bilateral(5);
  is_bilateral[0] = false;
  is_bilateral[1] = false;
  is_bilateral[2] = false;
  is_bilateral[3] = false;
  is_bilateral[4] = false;

  std::vector<float> b(5);
  b[0] =   2.0;
  b[1] =   2.0;
  b[2] =   2.0;
  b[3] =   2.0;
  b[4] =   2.0;

  std::vector<float> x(5);
  x[0] =  0.0;
  x[1] =  0.0;
  x[2] =  0.0;
  x[3] =  0.0;
  x[4] =  0.0;

  ProfilingData profiling_data;

  solve_eq_host(
                 5
                 , 13
                 , row
                 , column
                 , values
                 , x
                 , b
                 , params
                 , profiling_data
                 );

  std::cout << "x = [" << x[0];

  for(unsigned int i = 1u; i< 5u;++i)
    std::cout << " " << x[i];

  std::cout << " ] " << std::endl;

  x[0] =  0.0;
  x[1] =  0.0;
  x[2] =  0.0;
  x[3] =  0.0;
  x[4] =  0.0;

  solve_eq_device(
                   5
                   , 13
                   , row
                   , column
                   , values
                   , x
                   , b
                   , params
                   , profiling_data
                   );

  std::cout << "x = [" << x[0];

  for(unsigned int i = 1u; i< 5u;++i)
    std::cout << " " << x[i];
  
  std::cout << " ] " << std::endl;
}


inline void run_all_unit_tests()
{
//  test_minimum_map();
//  test_find_active();
//  test_project();
//  test_jacobian_operator();
//  test_jacobian_transpose_operator();
//  test_preconditioner_operator();
//  test_sparse_row();
//  test_psor();
//  test_schur_operator();
//  test_params();
//  test_min_map();
//  test_larger_min_map();
  test_solve_eq();
  test_solve_lcp();
}

// UNIT_TESTING_H
#endif
