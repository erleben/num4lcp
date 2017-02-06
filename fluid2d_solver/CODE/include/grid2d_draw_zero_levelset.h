#ifndef GRID2D_DRAW_ZERO_LEVELSET_H
#define GRID2D_DRAW_ZERO_LEVELSET_H

#include <grid2d.h>
#include <util_vec.h>
#include <util_svg_file.h>

namespace details
{

  template<typename P, typename T>
  inline void handle_triangle(Vec<2,P> const & p0
                                       , Vec<2,P> const & p1
                                       , Vec<2,P> const & p2
                                       , T const & d0
                                       , T const & d1
                                       , T const & d2
                                       , SVGFile * svg
                                       )
  {
    int const in0 = d0 < 0;  // inside/outside markers
    int const in1 = d1 < 0;
    int const in2 = d2 < 0;

    P const l0 = static_cast<P> ( - d0 / (d1-d0) );
    P const l1 = static_cast<P> ( - d1 / (d2-d1) );
    P const l2 = static_cast<P> ( - d2 / (d0-d2) );

    Vec<2,P> const e0 = p0 + l0*(p1-p0);
    Vec<2,P> const e1 = p1 + l1*(p2-p1);
    Vec<2,P> const e2 = p2 + l2*(p0-p2);

    int const mask = in0 | in1 << 1 | in2 << 2;

    switch (mask)
    {
      case 0: // 0 0 0 -- all outside
        break;
      case 7: // 1 1 1 -- all inside
        break;

      case 1:  // 0 0 1 -- p0 inside
        svg->draw_line(e0[0], e0[1], e2[0], e2[1]);
        break;

      case 2:  // 0 1 0 -- p1 inside
        svg->draw_line(e0[0], e0[1], e1[0], e1[1]);
        break;

      case 4:  // 1 0 0 -- p2 inside
        svg->draw_line(e1[0], e1[1], e2[0], e2[1]);
        break;

      case 3:  // 0 1 1 -- p0 and p1 inside
        svg->draw_line(e1[0], e1[1], e2[0], e2[1]);
        break;

      case 5:  // 1 0 1 -- p0 and p2 inside
        svg->draw_line(e0[0], e0[1], e1[0], e1[1]);
        break;

      case 6:  // 1 1 0 -- p1 and p2 inside
        svg->draw_line(e0[0], e0[1], e2[0], e2[1]);
        break;
    };
  }

}//namespace details


template<typename T, typename P>
inline void draw_zero_levelset( Grid2D<T,P> const & grid, SVGFile *svg = 0 )
{
  typedef Vec<2,P> V;

  int const J = grid.J()-1;
  int const I = grid.I()-1;

  for(int j = 0; j < J; ++j)
  {
    for(int i = 0; i < I; ++i)
    {
      T const f00 = grid(i  ,j  );
      T const f10 = grid(i+1,j  );
      T const f11 = grid(i+1,j+1);
      T const f01 = grid(i  ,j+1);

      V const p00 = grid.origin() + V(i*grid.dx(),j*grid.dx());
      V const p10 = p00           + V(  grid.dx(),          0);
      V const p11 = p00           + V(  grid.dx(),  grid.dx());
      V const p01 = p00           + V(          0,  grid.dx());

      int  const k    = (j*(I+1))+i;
      bool const flip = (k%2)==0;

      if(flip)
      {
        details::handle_triangle( p00, p10, p11, f00, f10, f11, svg );
        details::handle_triangle( p00, p11, p01, f00, f11, f01, svg );
      }
      else
      {
        details::handle_triangle( p01, p00, p10, f01, f00, f10, svg );
        details::handle_triangle( p01, p10, p11, f01, f10, f11, svg );
      }
    }
  }
  
}

// GRID2D_DRAW_ZERO_LEVELSET_H
#endif
