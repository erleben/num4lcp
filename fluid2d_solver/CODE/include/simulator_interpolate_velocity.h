#ifndef SIMULATOR_INTERPOLATE_VELOCITY_H
#define SIMULATOR_INTERPOLATE_VELOCITY_H

#include <grid2d_interpolate_value.h>

#include <simulator.h>

/**
 * Interpolate velocity from the MAC grid.
 */
inline typename Simulator::V interpolate_velocity(Simulator const & sim, typename Simulator::V const & position)
{
  // Interpolate the velocity from the u and v grids

  float const u_value = interpolate_value(sim.m_u, position);
  float const v_value = interpolate_value(sim.m_v, position );

  return Simulator::V(u_value, v_value);
}

// SIMULATOR_INTERPOLATE_VELOCITY_H
#endif