#ifndef GRID2D_COMPUTE_CELL_CENTER_POSITION_H
#define GRID2D_COMPUTE_CELL_CENTER_POSITION_H

#include <grid2d.h>
#include <util_vec.h>

/**
 *
 * @return   The position of the (i,j) cell center.
 */
template<typename T, typename P>
inline Vec<2,P> compute_cell_center_position(Grid2D<T,P> const & grid, int i, int j)
{
  typedef Vec<2,P> V;

  return V( i*grid.dx(), j*grid.dx() ) + grid.origin();
}

// GRID2D_COMPUTE_CELL_CENTER_POSITION_H
#endif
