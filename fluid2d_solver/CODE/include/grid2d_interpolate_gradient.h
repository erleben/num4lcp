#ifndef GRID2D_INTERPOLATE_GRADIENT_H
#define GRID2D_INTERPOLATE_GRADIENT_H

#include <grid2d.h>
#include <util_vec.h>

#include <cmath>
#include <cassert>

/**
 * Get value at any point inside the grid.
 */
template<typename T, typename P>
inline Vec<2,T> interpolate_gradient(Grid2D<T,P> const & grid, Vec<2,P> const & p)
{
  using std::floor;

  typedef Vec<2,T> Vt;
  typedef Vec<2,P> Vp;

  //--- Determine valid index range wrt. grid border
  int const i_min  = 0;
  int const i_max  = grid.I() - 2;
  int const j_min  = 0;
  int const j_max  = grid.J() - 2;

  //--- transform p into a local coordinate frame where grid origin is zero and dx is one.
  Vp const q =  (p-grid.origin())/grid.dx();

  //--- Determine the indices of four cell center containing p: (i,j), (i+1,j), (i,j+1), and (i+1,j+1)
  P const x = static_cast<P>(floor(q[0]));
  P const y = static_cast<P>(floor(q[1]));

  int i = (int)x;
  int j = (int)y;

  //--- Clamp to the borders of the grid
  P s = q[0] - x;
  if(i < i_min)
  {
    i = i_min;
    s = 0;
  }
  else if(i > i_max )
  {
    i = i_max;
    s = 1;
  }

  P t = q[1] - y;
  if(j < j_min)
  {
    j = j_min;
    t = 0;
  }
  else if(j > j_max )
  {
    j = j_max;
    t = 1;
  }

  //--- Now interpolate the gradient value at point p
  T const & v00 = grid(i  ,j  );
  T const & v10 = grid(i+1,j  );
  T const & v01 = grid(i  ,j+1);
  T const & v11 = grid(i+1,j+1);

  // Get gradient x-component by differentiating interpolation formula wrt. s,
  //
  //   ((1-t)*(s*v10 + (1-s)*v00)) + (t*(s*v11 + (1-s)*v01))
  //
  // and y-component by differentiation wrt. t. The result is
  return Vt(  (1-t)*(v10 - v00) + t*( v11 - v01)
           , s*(v11 - v10) + (1-s)*( v01 - v00)  );

}

// GRID2D_INTERPOLATE_GRADIENT_H
#endif
