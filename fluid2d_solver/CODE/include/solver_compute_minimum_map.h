#ifndef SOLVER_COMPUTE_MINIMUM_MAP_H
#define SOLVER_COMPUTE_MINIMUM_MAP_H

#include <cusp/array1d.h>
#include <cusp/copy.h>

#include <cassert>

namespace details
{

  struct MinimumMapFunctor
  {
    template <typename Tuple>
    __host__ __device__ void operator()(Tuple t)
    {
      //--- Tuple pattern: (y, x, is_bilateral, result )
      //---
      //--- if is_bilateral then
      //---
      //---   result = y
      //---
      //--- else
      //---
      //---   result = min(y,x)
      //---
      thrust::get<3>(t) = thrust::get<2>(t) ?  thrust::get<0>(t) : thrust::min(thrust::get<1>(t), thrust::get<0>(t));
    }
  };

  template <
  typename InputIterator1
  , typename InputIterator2
  , typename InputIterator3
  , typename OutputIterator
  >
  inline __host__ void compute_minimum_map(  InputIterator1 first_y
                                           , InputIterator1 last_y
                                           , InputIterator2 first_x
                                           , InputIterator3 first_is_bilateral
                                           , OutputIterator first_H
                                           )
  {
    size_t const N = last_y - first_y;

    thrust::for_each(
                     thrust::make_zip_iterator(thrust::make_tuple(first_y, first_x, first_is_bilateral, first_H)),
                     thrust::make_zip_iterator(thrust::make_tuple(first_y, first_x, first_is_bilateral, first_H)) + N,
                     MinimumMapFunctor()
                     );
  }

}// end namespace details


template <typename V, typename VB>
inline __host__ void compute_minimum_map(V const & y, V const & x, VB const & is_bilateral, V & H)
{
  assert(y.size() == x.size()            || !"compute_minimum_map(): x and y was not of same size");
  assert(x.size() == is_bilateral.size() || !"compute_minimum_map(): x and is_bilateral was not of same size");
  assert(x.size() == H.size()            || !"compute_minimum_map(): x and H was not of same size");
  
  details::compute_minimum_map(y.begin(), y.end(), x.begin(), is_bilateral.begin(), H.begin());
}


// SOLVER_COMPUTE_MINIMUM_MAP_H
#endif
