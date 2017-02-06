#ifndef GRID2D_EXTRAPOLATE_H
#define GRID2D_EXTRAPOLATE_H

#include <grid2d.h>

/**
 * Apply several iterations of a very simple "Jacobi"-style propagation
 * of valid velocity data in all directions.
 */
inline void extrapolate(Grid2D<float,float> & grid, Grid2D<char,float> & valid)
{
  Grid2D<char,float>  valid_old;
  Grid2D<float,float> grid_tmp;

  int const max_layers = 10;

  for(int layers = 0; layers < max_layers; ++layers)
  {
    valid_old = valid;
    grid_tmp  = grid;

    for(int j = 1; j < grid.J()-1; ++j)
      for(int i = 1; i < grid.I()-1; ++i)
      {
        float sum = 0;
        int count = 0;

        if(!valid_old(i,j))
        {

          if(valid_old(i+1,j))
          {
            sum += grid(i+1,j);
            ++count;
          }
          if(valid_old(i-1,j))
          {
            sum += grid(i-1,j);
            ++count;
          }
          if(valid_old(i,j+1))
          {
            sum += grid(i,j+1);
            ++count;
          }
          if(valid_old(i,j-1))
          {
            sum += grid(i,j-1);
            ++count;
          }

          //If any of neighbour cells were valid,
          //assign the cell their average value and tag it as valid
          if(count > 0)
          {
            grid_tmp(i,j) = sum /(float)count;
            valid(i,j) = 1;
          }

        }
      }
    grid = grid_tmp; //update with the new changes
  }  
}

// GRID2D_EXTRAPOLATE_H
#endif
