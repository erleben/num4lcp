#ifndef SOLVER_FIND_ACTIVE_H
#define SOLVER_FIND_ACTIVE_H

#include <cusp/array1d.h>
#include <cusp/verify.h>                        // Needed for cusp::assert_same_dimensions
#include <thrust/iterator/transform_iterator.h>


namespace details
{
  struct ActiveTest
  {
    template <typename Tuple>
    __host__ __device__ void operator()(Tuple t)
    {
      //--- Tuple pattern: (y, x, is_solid, result )
      //---
      //--- if is_solid then
      //---
      //---   result = true
      //---
      //--- else
      //---
      //---   result = y < x
      //---

      thrust::get<3>(t) = (thrust::get<0>(t) < thrust::get<1>(t)) || thrust::get<2>(t);
    }
  };


  template<
      typename InputIterator1
    , typename InputIterator2
    , typename InputIterator3
    , typename OutputIterator
    >
  inline __host__ void find_active(  InputIterator1 first_y
                                   , InputIterator1 last_y
                                   , InputIterator2 first_x
                                   , InputIterator3 first_is_solid
                                   , OutputIterator first_is_active
                                   )
  {
    size_t const N = last_y - first_y;

    thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(first_y, first_x, first_is_solid, first_is_active))
                     , thrust::make_zip_iterator(thrust::make_tuple(first_y, first_x, first_is_solid, first_is_active)) + N
                     , ActiveTest()
                     );
  }

}// end of namespace details


template <typename Vector1, typename Vector2>
inline __host__ void find_active(  Vector1 const & y
                                 , Vector1 const & x
                                 , Vector2 const & is_solid
                                 , Vector2 & is_active
                                 )
{
  cusp::assert_same_dimensions(y, x, is_solid, is_active);
  
  details::find_active(y.begin(), y.end(), x.begin(), is_solid.begin(), is_active.begin());
}


// SOLVER_FIND_ACTIVE_H
#endif
