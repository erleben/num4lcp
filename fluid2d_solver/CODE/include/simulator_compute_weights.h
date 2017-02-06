#ifndef SIMULATOR_COMPUTE_WEIGHTS_H
#define SIMULATOR_COMPUTE_WEIGHTS_H

#include <simulator.h>

#include <util_fraction_inside.h>


/**
 *
 * Compute finite-volume style face-weights for fluid from nodal signed distances
 *
 *
 * The weights are the fraction of the cell side not occupied by solid. Observe
 * that there could be both solid, air (vaccum) and liqiud fractions along a
 * cell side. The weights would correspond to both liqiud and air fractions.
 *
 */
inline void compute_weights(Simulator & sim)
{
  using std::min;
  using std::max;

  for(int j = 0; j < sim.m_u_weights.J(); ++j)
  {
    for(int i = 0; i < sim.m_u_weights.I(); ++i)
    {
      sim.m_u_weights(i,j) = 1 - fraction_inside(sim.m_phi_solid(i,j+1), sim.m_phi_solid(i,j));

      sim.m_u_weights(i,j) = max( min( 1.0f, sim.m_u_weights(i,j)), 0.0f);
    }
  }

  for(int j = 0; j < sim.m_v_weights.J(); ++j)
  {
    for(int i = 0; i < sim.m_v_weights.I(); ++i)
    {
      sim.m_v_weights(i,j) = 1 - fraction_inside(sim.m_phi_solid(i+1,j), sim.m_phi_solid(i,j));

      sim.m_v_weights(i,j) = max( min( 1.0f, sim.m_v_weights(i,j)), 0.0f);
    }
  }
}

// SIMULATOR_COMPUTE_WEIGHTS_H
#endif