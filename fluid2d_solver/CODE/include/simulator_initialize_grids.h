#ifndef SIMULATOR_INITIALIZE_GRIDS_H
#define SIMULATOR_INITIALIZE_GRIDS_H

#include <simulator.h>

/**
 *
 *
 * @param L   Characteristic length scale of the computational domain.
 * @param I   Number of cells along the x-axis
 * @param J   Number of cells along the y-axis
 */
inline void initialize_grids(Simulator & sim, float const & L, int const & I, int const & J)
{
  sim.m_I         = I;
  sim.m_J         = J;
  float const dx  = L / static_cast<float>(I-1);
  sim.m_dx        = dx;

  sim.m_particle_radius   = dx/sqrt(2.0f);
  sim.m_push_fraction     = 0.15f;

  sim.m_phi_liquid.initialize( I,   J,    0.0f,     0.0f,    dx, 0.0f);
  sim.m_phi_solid.initialize(  I+1, J+1, -dx/2.0f, -dx/2.0f, dx, 0.0f);

  sim.m_u.initialize(          I+1, J,   -dx/2.0f,  0.0f,    dx, 0.0f);
  sim.m_u_tmp.initialize(      I+1, J,   -dx/2.0f,  0.0f,    dx, 0.0f);
  sim.m_u_weights.initialize(  I+1, J,   -dx/2.0f,  0.0f,    dx, 0.0f);
  sim.m_u_valid.initialize(    I+1, J,   -dx/2.0f,  0.0f,    dx, 0   );

  sim.m_v.initialize(          I, J+1,    0.0f,    -dx/2.0f, dx, 0.0f);
  sim.m_v_tmp.initialize(      I, J+1,    0.0f,    -dx/2.0f, dx, 0.0f);
  sim.m_v_weights.initialize(  I, J+1,    0.0f,    -dx/2.0f, dx, 0.0f);
  sim.m_v_valid.initialize(    I, J+1,    0.0f,    -dx/2.0f, dx, 0   );
}

// SIMULATOR_INITIALIZE_GRIDS_H
#endif