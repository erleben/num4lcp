#ifndef SOLVER_FIND_FREE_H
#define SOLVER_FIND_FREE_H

#include <thrust/iterator/transform_iterator.h>


template<typename V>
inline __host__ void find_free(
                               V const & is_active
                               , V       & is_free
                               )
{
  typedef typename V::value_type T;

  thrust::transform(is_active.begin(), is_active.end(), is_free.begin(), thrust::logical_not<T>());
}


template<typename V>
inline __host__ const V find_free(V const & is_active)
{
  typedef typename V::value_type T;
  
  V is_free(is_active.size());

  find_free(is_active, is_free);

  return is_free;
}


// SOLVER_FIND_FREE_H
#endif
