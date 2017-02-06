#ifndef GRID2D_WRITE_TO_FILE_H
#define GRID2D_WRITE_TO_FILE_H

#include <grid2d.h>
#include <util_vec.h>

#include <iostream>
#include <fstream>
#include <sstream>

template<typename T, typename P>
inline void write_to_file( Grid2D<T,P> const & grid, std::string const & filename )
{
  std::ofstream file;
  file.open ( filename.c_str() );

  if(!file.is_open())
  {
    std::cerr << "could not open file " << filename << std::endl;
  }
  file << grid.I() << std::endl;
  file << grid.J() << std::endl;
  file << grid.dx() << std::endl;
  for (int j = 0; j < grid.J(); ++j)
  {
    for (int i = 0; i < grid.I(); ++i)
    {

      file << grid(i,j) << "\t";
    }
  }
  file << std::endl;
  file.flush();
  file.close();
}

// GRID2D_WRITE_TO_FILE_H
#endif
