#ifndef SIMULATOR_COMPUTE_CELL_HAS_LIQUID_H
#define SIMULATOR_COMPUTE_CELL_HAS_LIQUID_H

#include <simulator.h>

inline bool compute_cell_has_liquid(Simulator & sim, unsigned int const & i, unsigned int const & j)
{
//  unsigned int const ip1 = i+1u < sim.m_I ? i+1u :       i;
//  unsigned int const im1 = i    >      0u ? i-1u :      0u;
//  unsigned int const jp1 = j+1u < sim.m_J ? j+1u :       j;
//  unsigned int const jm1 = j    >      0u ? j-1u :      0u;
//
  float const phi_c   = sim.m_phi_liquid(i,j);
//  float const phi_r   = sim.m_phi_liquid(ip1,j);
//  float const phi_l   = sim.m_phi_liquid(im1,j);
//  float const phi_u   = sim.m_phi_liquid(i,jp1);
//  float const phi_d   = sim.m_phi_liquid(i,jm1);
//  float const phi_ur  = sim.m_phi_liquid(ip1,jp1);
//  float const phi_ul  = sim.m_phi_liquid(im1,jp1);
//  float const phi_dr  = sim.m_phi_liquid(ip1,jm1);
//  float const phi_dl  = sim.m_phi_liquid(im1,jm1);
//
//  float const psi_c   = phi_c;
//  float const psi_dl  = 0.25f*(phi_dl + phi_d  + phi_c  + phi_l);
//  float const psi_dr  = 0.25f*(phi_d  + phi_dr + phi_r  + phi_c);
//  float const psi_ur  = 0.25f*(phi_c  + phi_r  + phi_ur + phi_u);
//  float const psi_ul  = 0.25f*(phi_l  + phi_c  + phi_u  + phi_ul);
//
  //return ( psi_dl < 0 || psi_dr < 0 || psi_ur < 0 || psi_ul < 0 || psi_c < 0);  // makes velocties funny?
  return (  phi_c < 0);
}

// SIMULATOR_COMPUTE_CELL_HAS_LIQUID_H
#endif
