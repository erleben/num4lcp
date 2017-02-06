#ifndef GRID2D_H
#define GRID2D_H

#include <util_vec.h>

#include <vector>
#include <cassert>

/**
 * Data structure for a uniform regular grid.
 *
 *   (I,J) = (5,4) => ie. five cells wide and four cells tall
 *
 *          +-----+-----+-----+-----+-----+
 *          |     |     |     |     |     |
 *          |  *  |  *  |  *  |  *  |  *  |
 *          |(0,3)|(1,3)|(2,3)|(3,3)|(4,3)|
 *          +-----+-----+-----+-----+-----+
 *          |     |     |     |     |     |
 *          |  *  |  *  |  *  |  *  |  *  |
 *          |(0,2)|(1,2)|(2,2)|(3,2)|(4,2)|
 *          +-----+-----+-----+-----+-----+
 *          |     |     |     |     |     |
 *          |  *  |  *  |  *  |  *  |  *  |
 *          |(0,1)|(1,1)|(2,1)|(3,1)|(4,1)|
 *          +-----+-----+-----+-----+-----+
 *          |     |     |     |     |     |
 * origin --|->*  |  *  |  *  |  *  |  *  |  dy = dx
 *          |(0,0)|(1,0)|(2,0)|(3,0)|(4,0)|
 *          +-----+-----+-----+-----+-----+
 *                                     dx
 *
 * @tparam T   The data type stored in the grid cell centers. 
 * @tparam P   The precision used for positions, sizes etc.
 */
template<typename T, typename P>
class Grid2D
{
public:

  template<typename T2, typename P2>
  friend T2 interpolate_value(Grid2D<T2,P2> const & grid, Vec<2,P2> const & p);

  template<typename T2, typename P2>
  friend Vec<2,T2> interpolate_gradient(Grid2D<T2,P2> const & grid, Vec<2,P2> const & p);

  template<typename T2, typename P2>
  friend Vec<2,P2> compute_cell_center_position(Grid2D<T2,P2> const & grid,int i, int j);

  template<typename T2, typename P2>
  friend void compute_cell_indices(Grid2D<T2,P2> const & grid, Vec<2,P2> const & p, int & i, int & j);

  template<typename T2, typename P2>
  friend void draw_zero_levelset( Grid2D<T2,P2> const & grid );

public:

  typedef Vec<2,P> V;

  typedef typename std::vector<T>::iterator       iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

public:

  int            m_I;      ///< Number of cells along x-axis.
  int            m_J;      ///< Number of cells along y-axis.
  V              m_origin; ///< Origin of grid, lower-left cell center.
  P              m_dx;     ///< The cell size (length of cell side).
  std::vector<T> m_data;   ///< The data values of the cell centers stored as one linear array (row-by-row manner).

public:


  int I()      const { return this->m_I;      }
  int J()      const { return this->m_J;      }
  P   dx()     const { return this->m_dx;     }
  V   origin() const { return this->m_origin; }

  iterator       begin()       { return this->m_data.begin(); }
  const_iterator begin() const { return this->m_data.begin(); }
  iterator       end()         { return this->m_data.end();   }
  const_iterator end()   const { return this->m_data.end();   }

public:

  Grid2D()
  : m_I(0)
  , m_J(0)
  , m_origin(0,0)
  , m_dx(0)
  , m_data()
  {}

public:

  /**
   * Initialize grid.
   *
   * @param I      The number of cells along the x-axis.
   * @param J      The number of cells along the y-axis.
   * @param x      The x-coordinate of the lower-left cell center of the grid.
   * @param y      The y-coordinate of the lower-left cell center of the grid.
   * @param dx     The horizontal/vertical spacing between cell centers (same as cell length/size).
   * @param value  The default value to fill into the cell centers.
   */
  void initialize(int I, int J, P const & x, P const & y, P const & dx, T const & value)
  {
    assert(I>=0 || !"Illegal I value");
    assert(J>=0 || !"Illegal J value");
    assert(dx>0 || !"Illegal dx value");

    this->m_I      = I;
    this->m_J      = J;
    this->m_origin = V(x,y);
    this->m_dx     = dx;
    this->m_data.resize(I*J,value);
  }

  size_t size() const
  {
    return this->m_data.size();
  }

  T const & operator()(int i, int j) const
  {
    assert( this->m_I > 0 || !"Illegal usage, grid is empty");
    assert( this->m_J > 0 || !"Illegal usage, grid is empty");
    assert(i >= 0         || !"illegal i value");
    assert(i < this->m_I  || !"illegal i value");
    assert(j >= 0         || !"illegal j value");
    assert(j < this->m_J  || !"illegal j value");

    return this->m_data[i + this->m_I*j];
  }

  T & operator()(int i, int j)
  {
    assert( this->m_I > 0 || !"Illegal usage, grid is empty");
    assert( this->m_J > 0 || !"Illegal usage, grid is empty");
    assert(i >= 0         || !"illegal i value");
    assert(i < this->m_I  || !"illegal i value");
    assert(j >= 0         || !"illegal j value");
    assert(j < this->m_J  || !"illegal j value");

    return this->m_data[i + this->m_I*j];
  }

};

// GRID2D_H
#endif
