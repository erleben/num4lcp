#ifndef SOLVER_SPARSE_ROW_DOT_H
#define SOLVER_SPARSE_ROW_DOT_H

#include <cusp/csr_matrix.h>
#include <cusp/array1d.h>
#include <cusp/blas/blas.h>


template <typename matrix_type, typename vector_type>
inline __host__ void sparse_row_dot(
                             matrix_type                     const & A
                             , vector_type                     const & x
                             , typename vector_type::size_type const & row_index
                             , typename vector_type::value_type      & dot_product
                             )
{
  typedef typename vector_type::value_type real_type;
  
  dot_product = real_type(0);
  
  size_t const num_elem_row = A.row_offsets[row_index+1] - A.row_offsets[row_index];
  if (num_elem_row == 0)
  {
    return;
  }
  
  size_t idx;
  size_t offset;
  for (size_t i = 0; i < num_elem_row; i++)
  {
    offset = A.row_offsets[row_index] + i;
    idx = A.column_indices[offset];
    dot_product += A.values[offset] * x[idx];
  }
  
}


// SOLVER_SPARSE_ROW_DOT_H
#endif
