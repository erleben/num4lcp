#ifndef SIMULATOR_RUN_H
#define SIMULATOR_RUN_H

#include <grid2d_extrapolate.h>

#include <simulator_compute_cfl.h>
#include <simulator_compute_volume.h>
#include <simulator_compute_energy.h>
#include <simulator_compute_momentum.h>
#include <simulator_add_gravity_force.h>
#include <simulator_apply_boundary_constraint_on_velocity.h>
#include <simulator_advect_velocity.h>
#include <simulator_advect_particles.h>
#include <simulator_compute_phi_liquid.h>
#include <simulator_compute_weights.h>
#include <simulator_solve_pressure_poisson_equation.h>



inline void run(Simulator & sim, Params const & params)
{
  float const T = params.time_step();

  float t = 0;

  while(t < T)
  {
    START_TIMER("total");

    //--- Compute how big a time step we can take -----------------------------
    float dt = compute_cfl( sim );

    if (t + dt > T)
    {
      dt = T - t;
    }
    
    //--- Deal with advection -------------------------------------------------
    START_TIMER("advect_particles");
    advect_particles( sim, dt);
    STOP_TIMER("advect_particles");

    START_TIMER("advect_velocity");
    advect_velocity( sim, dt);
    STOP_TIMER("advect_velocity");

    //--- Integrate gravity force ---------------------------------------------
    START_TIMER("gravity");
    add_gravity_force( sim, dt);
    STOP_TIMER("gravity");

    //--- Perform pressure projection -----------------------------------------
    //Estimate the liquid signed distance
    START_TIMER("compute_phi");
    compute_phi_liquid( sim );
    STOP_TIMER("compute_phi");

    // Compute finite-volume type face area weight for each velocity sample.
    START_TIMER("compute_weights");
    compute_weights(sim);
    STOP_TIMER("compute_weights");

    // Set up and solve the variational pressure solve.
    solve_pressure_poisson_equation(sim, dt, params);

    // Pressure projection only produces valid velocities in faces with non-zero
    // associated face area. Because the advection step may interpolate from
    // these invalid faces, we must extrapolate velocities from the fluid domain
    // into these zero-area faces.
    START_TIMER("extrapolate");
    extrapolate(sim.m_u, sim.m_u_valid);
    extrapolate(sim.m_v, sim.m_v_valid);
    STOP_TIMER("extrapolate");

    // For extrapolated velocities, replace the normal component with
    // that of the object.
    START_TIMER("apply_boundary_constraint");
    apply_boundary_constraint_on_velocity( sim );
    STOP_TIMER("apply_boundary_constraint");

    //--- Advance time to next time step ------------------------------------
    t += dt;

    STOP_TIMER("total");

    if(params.profiling())
    {
      float       volume_error   = 0.0f;
      float const volume         = compute_volume(sim, volume_error);
      float       energy_error   = 0.0f;
      float const energy         = compute_energy(sim, energy_error);
      float       momentum_error = 0.0f;
      float const momentum       = compute_momentum(sim, momentum_error);

      RECORD("timestep", dt);
      RECORD("energy", energy);
      RECORD("volume", volume);
      RECORD("momentum", momentum);
      RECORD("energy_error", energy_error);
      RECORD("volume_error", volume_error);
      RECORD("momentum_error", momentum_error);
    }

  }
}

// SIMULATOR_RUN_H
#endif
