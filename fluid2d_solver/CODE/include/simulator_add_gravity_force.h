#ifndef SIMULATOR_ADD_GRAVITY_FORCE_H
#define SIMULATOR_ADD_GRAVITY_FORCE_H

#include <simulator.h>

inline void add_gravity_force(Simulator & sim, float const & dt)
{
  for(int j = 0; j < sim.m_v.J(); ++j)
  {
    for(int i = 0; i < sim.m_v.I(); ++i)
    {
      // Kenny:
      //
      //   F dt dV / rho = F dt / dm = dm G dt /dm = G dt
      //
      // Because rho = dm/dV, here dm is the liquid mass.

      sim.m_v(i,j) -= 9.82f * dt;
    }
  }
}

// SIMULATOR_ADD_GRAVITY_FORCE_H
#endif