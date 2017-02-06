#ifndef SIMULATOR_COMPUTE_MOMENTUM_H
#define SIMULATOR_COMPUTE_MOMENTUM_H

#include <simulator.h>
#include <simulator_interpolate_velocity.h>
#include <grid2d_compute_cell_center_position.h>

/**
 *
 * Compute approximate magnitude of total linear momentum  of liquid.
 */
inline float compute_momentum(Simulator & sim, float & error)
{
  typedef typename Simulator::V V;

  float const  dv    = sim.m_dx*sim.m_dx;

  V M     = V(0.0f, 0.0f); // Total linear moementum
    error = 0.0f;


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
      V     const dM          = dv*interpolate_velocity( sim, p );

      M = (phi_center < 0)? M + dM : M;

      error = (phi_center < 0 && (phi_left>0 || phi_right > 0 || phi_up > 0 || phi_down > 0)) ? error + mag(dM) : error;
    }
  }

  return mag(M);
}

// SIMULATOR_COMPUTE_MOMENTUM_H
#endif