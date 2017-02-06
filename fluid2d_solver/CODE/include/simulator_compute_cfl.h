#ifndef SIMULATOR_COMPUTE_CFL_H
#define SIMULATOR_COMPUTE_CFL_H

#include <simulator.h>

#include <cmath>


inline float compute_cfl(Simulator const & sim)
{
  typedef typename Simulator::Grid2Df::const_iterator const_iterator;

  using std::max;
  using std::fabs;

  float max_vel = 0;

  for(const_iterator iter = sim.m_u.begin(); iter != sim.m_u.end(); ++iter)
    max_vel = max(max_vel, fabs(*iter));

  for( const_iterator iter = sim.m_v.begin(); iter != sim.m_v.end(); ++iter)
    max_vel = max(max_vel, fabs(*iter));

  return sim.m_dx / max_vel;
}

// SIMULATOR_COMPUTE_CFL_H
#endif