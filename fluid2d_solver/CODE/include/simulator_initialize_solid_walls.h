#ifndef SIMULATOR_INITIALIZE_SOLID_WALLS_H
#define SIMULATOR_INITIALIZE_SOLID_WALLS_H

#include <simulator.h>

/**
 * Initialize the grid-based signed distance field that
 * dictates the position of the solid wall boundary.
 *
 * @param phi    A signed distance field expression that represent the
 *               solid geometry and can be used to evaluate the signed
 *               distance field on the solid distance field grid nodes.
 */
template<typename F>
inline void initialize_solid_walls(Simulator & sim, DistanceFieldExp<F> const & phi)
{
  for(int j = 0; j < sim.m_phi_solid.J(); ++j)
  {
    for(int i = 0; i < sim.m_phi_solid.I(); ++i)
    {
      typename Simulator::V const pos = compute_cell_center_position(sim.m_phi_solid,i,j);
      sim.m_phi_solid(i,j) = phi(pos);
    }
  }
}

// SIMULATOR_INITIALIZE_SOLID_WALLS_H
#endif