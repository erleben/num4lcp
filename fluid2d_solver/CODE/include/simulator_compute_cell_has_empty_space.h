#ifndef SIMULATOR_COMPUTE_CELL_HAS_EMPTY_SPACE_H
#define SIMULATOR_COMPUTE_CELL_HAS_EMPTY_SPACE_H

#include <simulator.h>

inline bool compute_cell_has_empty_space(Simulator & sim, unsigned int const & i, unsigned int const & j)
{
  float const psi_dl  = sim.m_phi_solid(  i,  j);
  float const psi_dr  = sim.m_phi_solid(i+1,  j);
  float const psi_ul  = sim.m_phi_solid(  i,j+1);
  float const psi_ur  = sim.m_phi_solid(i+1,j+1);

  return ( psi_dl > 0 || psi_dr > 0 || psi_ur > 0 || psi_ul > 0 );
}

// SIMULATOR_COMPUTE_CELL_HAS_EMPTY_SPACE_H
#endif
