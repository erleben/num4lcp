#ifndef SIMULATOR_COMPUTE_PHI_LIQUID_H
#define SIMULATOR_COMPUTE_PHI_LIQUID_H

#include <grid2d_compute_cell_indices.h>
#include <grid2d_compute_cell_center_position.h>

#include <simulator.h>

inline void compute_phi_liquid(Simulator & sim)
{
  using std::max;
  using std::min;

  typedef typename Simulator::V V;

  //--- Estimate signed distance field from particles
  float const init_value = 3*sim.m_dx;
  std::fill( sim.m_phi_liquid.begin(), sim.m_phi_liquid.end(), init_value );

  float const R = 1.02f * sim.m_particle_radius;
  for(size_t p = 0; p < sim.m_particles.size(); ++p)
  {
    V const point = sim.m_particles[p];

    int i;
    int j;
    compute_cell_indices(sim.m_phi_liquid, point, i, j);

    //--- compute distance to surrounding few points, ---------------------
    //--- keep if it's the minimum ----------------------------------------
    int   const i_min = max(i-2,0);
    int   const j_min = max(j-2,0);
    int   const i_max = min(i+2,sim.m_phi_liquid.I()-1);
    int   const j_max = min(j+2,sim.m_phi_liquid.J()-1);

    for(int j_off = j_min; j_off <= j_max; ++j_off)
      for(int i_off = i_min; i_off <= i_max; ++i_off)
      {
        V const pos = compute_cell_center_position(sim.m_phi_liquid, i_off, j_off);

        float const phi_guess = dist(pos, point) - R;

        sim.m_phi_liquid(i_off,j_off) = min(sim.m_phi_liquid(i_off,j_off), phi_guess);
      }
  }

}

// SIMULATOR_COMPUTE_PHI_LIQUID_H
#endif