#ifndef SOLVER_PROJECT_H
#define SOLVER_PROJECT_H


#include <cusp/array1d.h>
#include <cusp/copy.h>

#include <cassert>


namespace details
{

  template <typename T>
  struct MaxFunctor
  {
    template <typename Tuple>
    __host__ __device__ void operator()(Tuple t)
    {
      //--- Tuple pattern: (x, is_bilateral, result )
      //---
      //--- if is_bilateral then
      //---
      //---   result = x
      //---
      //--- else
      //---
      //---   result = max(0,x)
      //---
      thrust::get<2>(t) = thrust::get<1>(t)  ? thrust::get<0>(t)  : thrust::max(T(0), thrust::get<0>(t));
    }
  };


  template <
  typename InputIterator1
  , typename InputIterator2
  , typename OutputIterator
  >
  inline __host__ void project(  InputIterator1 first_x
                               , InputIterator1 last_x
                               , InputIterator2 first_is_bilateral
                               , OutputIterator first_project_x
                               )
  {
    typedef typename InputIterator1::value_type T;

    size_t const N = last_x - first_x;

    thrust::for_each(
                     thrust::make_zip_iterator(thrust::make_tuple(first_x, first_is_bilateral, first_project_x))
                     , thrust::make_zip_iterator(thrust::make_tuple(first_x, first_is_bilateral, first_project_x)) + N
                     , MaxFunctor<T>()
                     );
  }

}// end of namespace details


template <typename V, typename VB>
inline __host__ void project(  V const & x, VB const & is_bilateral, V & projected_x)
{
  assert(x.size() == projected_x.size()  || !"project(): x and projected_x incompatible dimensions");
  assert(x.size() == is_bilateral.size() || !"project(): x and is_bilateral incompatible dimensions");

  details::project(
                   x.begin()
                   , x.end()
                   , is_bilateral.begin()
                   , projected_x.begin()
                   );
}


// SOLVER_PROJECT_H
#endif
