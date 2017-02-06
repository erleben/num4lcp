#ifndef GRID2D_COMPUTE_CELL_INDICES_H
#define GRID2D_COMPUTE_CELL_INDICES_H

#include <grid2d.h>
#include <util_vec.h>

#include <cmath>
#include <cassert>

/**
 * Get value at any point inside the grid.
 */
template<typename T, typename P>
inline void compute_cell_indices(Grid2D<T,P> const & grid
                                 , Vec<2,P> const & p
                                 , int & i
                                 , int & j
                                 )
{
  using std::min;
  using std::max;
  using std::floor;

  typedef Vec<2,P> V;

  //--- Determine valid index range wrt. grid border
  int const i_min  = 0;
  int const i_max  = grid.I() - 2;
  int const j_min  = 0;
  int const j_max  = grid.J() - 2;

  //--- transform p into a local coordinate frame where grid origin is zero and dx is one.
  V const q =  (p-grid.origin())/grid.dx();

  //--- Determine the indices of four cell centers containing p: (i,j), (i+1,j), (i,j+1), and (i+1,j+1)
  P const x = static_cast<P>(floor(q[0]));
  P const y = static_cast<P>(floor(q[1]));

  i = (int)x;
  j = (int)y;

  i = max( min( i, i_max ), i_min );
  j = max( min( j, j_max ), j_min );
}

// GRID2D_COMPUTE_CELL_INDICES_H
#endif
