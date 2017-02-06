#ifndef SIMULATOR_ADVECT_VELOCITY_H
#define SIMULATOR_ADVECT_VELOCITY_H

#include <grid2d_compute_cell_center_position.h>

#include <simulator.h>
#include <simulator_compute_rk2.h>
#include <simulator_interpolate_velocity.h>

/**
 * Basic first order semi-Lagrangian advection of velocities
 */
inline void advect_velocity(Simulator & sim, float dt)
{
  typedef typename Simulator::V V;

  //--- semi-Lagrangian advection on u-component of velocity --------
  for(int j = 0; j < sim.m_u.J(); ++j)
    for(int i = 0; i < sim.m_u.I(); ++i)
    {
      V pos            = compute_cell_center_position(sim.m_u,i,j);
      pos              = compute_rk2( sim, pos, -dt);
      sim.m_u_tmp(i,j) = interpolate_velocity( sim, pos)[0];
    }

  //--- semi-Lagrangian advection on v-component of velocity -------
  for(int j = 0; j < sim.m_v.J(); ++j)
    for(int i = 0; i < sim.m_v.I(); ++i)
    {
      V pos            = compute_cell_center_position(sim.m_v,i,j);
      pos              = compute_rk2( sim, pos, -dt);
      sim.m_v_tmp(i,j) = interpolate_velocity( sim, pos)[1];
    }


  //--- move update velocities into u/v vectors ---------------------
  sim.m_u = sim.m_u_tmp;
  sim.m_v = sim.m_v_tmp;
}

// SIMULATOR_ADVECT_VELOCITY_H
#endif