#ifndef SIMULATOR_ADVECT_PARTICLES_H
#define SIMULATOR_ADVECT_PARTICLES_H

#include <grid2d_interpolate_value.h>
#include <grid2d_interpolate_gradient.h>

#include <simulator.h>
#include <simulator_compute_rk2.h>
#include <simulator_interpolate_velocity.h>

/**
 * Perform 2nd order Runge Kutta to move the particles in the fluid
 */
inline void advect_particles(Simulator & sim, float dt)
{
  typedef typename Simulator::V V;

  for(size_t p = 0; p < sim.m_particles.size(); ++p)
  {
    V const before     = sim.m_particles[p];
    sim.m_particles[p] = compute_rk2( sim, before, dt );
    V const after      = sim.m_particles[p];

    if(dist(before,after) > 3*sim.m_dx)
    {
      std::cerr << "Warning potential instability -- step size likely too large" << std::endl;
    }

    //--- Particles can still occasionally leave the domain due to truncation
    //--- errors, interpolation error, or large timesteps, so we project them
    //--- back in for good measure.

    float const phi_value = interpolate_value(sim.m_phi_solid, sim.m_particles[p]);
    if(phi_value < sim.m_particle_radius*sim.m_push_fraction)
    {
      V const normal = normalized(interpolate_gradient(sim.m_phi_solid, sim.m_particles[p]));
      sim.m_particles[p] -= (phi_value-sim.m_particle_radius)*normal;
    }
  }
}

// SIMULATOR_ADVECT_PARTICLES_H
#endif