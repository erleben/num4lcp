#ifndef SIMULATOR_COMPUTE_VOLUME_H
#define SIMULATOR_COMPUTE_VOLUME_H

#include <simulator.h>

/**
 *
 * Compute approximate volume of liqued.
 */
inline float compute_volume(Simulator & sim, float & error)
{
  float const dv = sim.m_dx*sim.m_dx;

  unsigned int N = 0u; // Number of cells close to the surface;
  float V = 0.0f;      // Total liquid volume sofar

  for(int j = 1; j < sim.m_J-1; ++j)
  {
    for(int i = 1; i < sim.m_I-1; ++i)
    {
      float const phi_center  = sim.m_phi_liquid(i,j);
      float const phi_right   = sim.m_phi_liquid(i+1,j);
      float const phi_left    = sim.m_phi_liquid(i-1,j);
      float const phi_up      = sim.m_phi_liquid(i,j+1);
      float const phi_down    = sim.m_phi_liquid(i,j-1);

      V = (phi_center < 0)? V+dv : V;
      N = (phi_center < 0 && (phi_left>0 || phi_right > 0 || phi_up > 0 || phi_down > 0)) ? N+1u : N;
    }
  }

  error = dv*N; // Conservative upper bound on how much volume could have been "miscalculated"
  return V;
}

// SIMULATOR_COMPUTE_VOLUME_H
#endif