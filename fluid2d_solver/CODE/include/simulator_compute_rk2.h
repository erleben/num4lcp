#ifndef SIMULATOR_COMPUTE_RK2_H
#define SIMULATOR_COMPUTE_RK2_H

#include <simulator.h>

/**
 * Apply RK2 to advect a point in the domain.
 */
inline typename Simulator::V compute_rk2(Simulator const & sim
                                  , typename Simulator::V const & position
                                  , float const & dt)
{
  typedef typename Simulator::V V;

  V input    = position;
  V velocity = interpolate_velocity( sim, input );
  velocity   = interpolate_velocity( sim, input + 0.5f*dt*velocity );
  input     += dt*velocity;
  return input;
}

// SIMULATOR_COMPUTE_RK2_H
#endif