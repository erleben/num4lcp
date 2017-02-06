#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <grid2d.h>
#include <util_distance_field.h>
#include <util_vec.h>

#include <vector>

class Simulator
{
public:

  typedef Vec<2,float> V;

  typedef Grid2D<float, float> Grid2Df;
  typedef Grid2D<char, float>  Grid2Dc;

public:

  int   m_I;  ///< Number of cells along the x-axis
  int   m_J;  ///< Number of cells along the y-axis
  float m_dx; ///< The length of a cell face (side length). Assumed to be uniform for both x and y directions.

  float            m_particle_radius;    ///< Radius of marker particles.
  float            m_push_fraction;      ///< Amount of particle radius that liquid is allowed to overlap solid walls.
  std::vector< V > m_particles;          ///< Storage for marker particle simulation.

  Grid2Df m_u;          ///< x-component of velocity field stored on vertical face centers
  Grid2Df m_v;          ///< y-component of velocity field stored on horizontal face centers
  Grid2Df m_u_tmp;      ///< temporary x-component of velocity field, used for flip-flopping
  Grid2Df m_v_tmp;      ///< temporary y-component of velocity field, used for flip-flopping
  Grid2Df m_u_weights;  ///< ??
  Grid2Df m_v_weights;  ///< ??
  Grid2Dc m_u_valid;    ///< ??
  Grid2Dc m_v_valid;    ///< ??

  Grid2Df m_phi_solid;  ///< The signed distance map of the solid wall boundaries.
  Grid2Df m_phi_liquid; ///< The signed distance map of the surface of the liquid.

public:

  unsigned int m_N;                               ///< Number of variables
  unsigned int m_K;                               ///< Number of non zeros in coefficient matrix
  std::vector<unsigned int> m_row_indices;        ///< Row indices of values in coefficient matrix
  std::vector<unsigned int> m_column_indices;     ///< Column indices of values in coefficient matrix
  std::vector<float>        m_values;             ///< Values in coefficient matrix
  std::vector<float>        m_x;                  ///< The solution vector
  std::vector<float>        m_b;                  ///< The right hand side vector
  std::vector<bool>         m_is_bilateral;       ///< Forced bilateral contraints

};

// SIMULATOR_H
#endif
