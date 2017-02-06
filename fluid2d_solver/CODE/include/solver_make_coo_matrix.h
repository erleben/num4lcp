#ifndef SOLVER_MAKE_COO_MATRIX_H
#define SOLVER_MAKE_COO_MATRIX_H

#include <gpu_check_error.h>

#include <thrust/copy.h>
#include <cusp/coo_matrix.h>

#include <cassert>
#include <vector>

template<typename I, typename T, typename M>
inline __host__  cusp::coo_matrix<I,  T,  M>  make_coo_matrix(
                                                              unsigned int const & n
                                                              , unsigned int const & k
                                                              , std::vector<unsigned int> const & row_indices
                                                              , std::vector<unsigned int> const & column_indices
                                                              , std::vector<float       > const & values
                                                              )
{
  CHECK_ERROR("make coo matrix invoked");

  assert(row_indices.size() >= k    || !"make_coo_matrix(): Internal error row indices did not have proper size");
  assert(column_indices.size() >= k || !"make_coo_matrix(): Internal error column indices did not have proper size");
  assert(values.size() >= k         || !"make_coo_matrix(): Internal error values did not have proper size");

  cusp::coo_matrix<I,  T,  M>  Acoo(n, n, k);
  CHECK_ERROR("Allocated space for Acoo matrix");

  thrust::copy(row_indices.begin(),    row_indices.begin()+k,    Acoo.row_indices.begin()    );
  thrust::copy(column_indices.begin(), column_indices.begin()+k, Acoo.column_indices.begin() );
  thrust::copy(values.begin(),         values.begin()+k,         Acoo.values.begin()         );
  CHECK_ERROR("Done copying data intoAcoo matrix");

  Acoo.sort_by_row_and_column();
  CHECK_ERROR("Sorted data of Acoo matrix");

  return Acoo;
}

// SOLVER_MAKE_COO_MATRIX_H
#endif
