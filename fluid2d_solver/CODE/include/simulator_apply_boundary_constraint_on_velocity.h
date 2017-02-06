#ifndef SIMULATOR_APPLY_BOUNDARY_CONSTRAINT_ON_VELOCITY_H
#define SIMULATOR_APPLY_BOUNDARY_CONSTRAINT_ON_VELOCITY_H

#include <grid2d_interpolate_gradient.h>
#include <grid2d_compute_cell_center_position.h>

#include <simulator.h>
#include <simulator_interpolate_velocity.h>

/**
 * For extrapolated points, replace the normal component
 * of velocity with the object velocity (in this case zero).
 */
inline void apply_boundary_constraint_on_velocity(Simulator & sim)
{
  typedef typename Simulator::V V;
  
  sim.m_u_tmp = sim.m_u;
  sim.m_v_tmp = sim.m_v;

  //
  // At lower grid resolutions, the normal estimate from the signed
  // distance function is poor, so it doesn't work quite as well.
  // An exact normal would do better.
  //

  //--- Constrain u -----------------
  for(int j = 0; j < sim.m_u.J(); ++j)
    for(int i = 0; i < sim.m_u.I(); ++i)
    {
      if(sim.m_u_weights(i,j) == 0) // Full solid wall horizontal edge
      {
        //apply boundary constraint
        V     const pos            = compute_cell_center_position(sim.m_u, i, j);
        V           vel            = interpolate_velocity( sim, pos);
        V     const normal         = normalized( interpolate_gradient( sim.m_phi_solid, pos ) );
        float const perp_component = dot(vel, normal);

        vel -= perp_component*normal;

        sim.m_u_tmp(i,j) = vel[0];
      }
    }

  //--- Constrain v ------------------
  for(int j = 0; j < sim.m_v.J(); ++j)
    for(int i = 0; i < sim.m_v.I(); ++i)
    {
      if(sim.m_v_weights(i,j) == 0)  // Full solid wall vertical edge
      {
        //apply boundary constraint
        V     const pos            = compute_cell_center_position(sim.m_v, i, j);
        V           vel            = interpolate_velocity( sim, pos);
        V     const normal         = normalized( interpolate_gradient( sim.m_phi_solid, pos ) );
        float const perp_component = dot(vel, normal);

        vel -= perp_component*normal;

        sim.m_v_tmp(i,j) = vel[1];
      }
    }

  //--- Update ----------------------
  sim.m_u = sim.m_u_tmp;
  sim.m_v = sim.m_v_tmp;
}

// SIMULATOR_APPLY_BOUNDARY_CONSTRAINT_ON_VELOCITY_H
#endif