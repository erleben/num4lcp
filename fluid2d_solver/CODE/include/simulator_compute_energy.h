#ifndef SIMULATOR_COMPUTE_ENERGY_H
#define SIMULATOR_COMPUTE_ENERGY_H

#include <simulator.h>
#include <simulator_interpolate_velocity.h>
#include <grid2d_compute_cell_center_position.h>

/**
 *
 * Compute approximate kinetic energy of liquid.
 */
inline float compute_energy(Simulator & sim, float & error)
{
  typedef typename Simulator::V V;

  float const  dv    = sim.m_dx*sim.m_dx;
  float        K     = 0.0f;      // Total kinetic energy of liquid
               error = 0.0f;      // Conservative upper bound on how much kinetic energy could have been "miscalculated"

  for(int j = 1; j < sim.m_J-1; ++j)
  {
    for(int i = 1; i < sim.m_I-1; ++i)
    {
      float const phi_center  = sim.m_phi_liquid(i,j);
      float const phi_right   = sim.m_phi_liquid(i+1,j);
      float const phi_left    = sim.m_phi_liquid(i-1,j);
      float const phi_up      = sim.m_phi_liquid(i,j+1);
      float const phi_down    = sim.m_phi_liquid(i,j-1);
      V     const p           = compute_cell_center_position(sim.m_phi_liquid,i,j);
      V     const u           = interpolate_velocity( sim, p );
      float const dK          = dv*dot(u,u)/2.0f;

      K = (phi_center < 0)? K+ dK : K;

      error = (phi_center < 0 && (phi_left>0 || phi_right > 0 || phi_up > 0 || phi_down > 0)) ? error + dK : error;
    }
  }

  return K;
}

// SIMULATOR_COMPUTE_ENERGY_H
#endif