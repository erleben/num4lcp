#ifndef GRID2D_DRAW_GRID_LINES_H
#define GRID2D_DRAW_GRID_LINES_H

#include <grid2d.h>
#include <util_vec.h>
#include <util_svg_file.h>


template<typename T, typename P>
inline void draw_grid_lines( Grid2D<T,P> const & grid, SVGFile * svg )
{
  typedef Vec<2,P> V;

  P const width  = (grid.I()-1u)*grid.dx();
  P const height = (grid.J()-1u)*grid.dx();

  for(int i = 0; i < grid.I(); ++i)
  {
    V const a(i*grid.dx(), 0     );
    V const b(i*grid.dx(), height);
    V const from = (grid.origin()+a);
    V const to   = (grid.origin()+b);

    svg->draw_line(from[0], from[1], to[0], to[1] );
  }

  for(int j = 0; j < grid.J(); ++j)
  {
    V const a(    0, j*grid.dx());
    V const b(width, j*grid.dx());
    V const from = (grid.origin()+a);
    V const to   = (grid.origin()+b);

    svg->draw_line(from[0], from[1], to[0], to[1] );
  }
}

// GRID2D_DRAW_GRID_LINES_H
#endif
