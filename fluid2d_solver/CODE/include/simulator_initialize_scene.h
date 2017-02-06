#ifndef SIMULATOR_INITIALIZE_SCENE_H
#define SIMULATOR_INITIALIZE_SCENE_H

#include <climits>

#include <simulator.h>
#include <simulator_initialize_grids.h>
#include <simulator_initialize_solid_walls.h>

#include <util_distance_field.h>


namespace details
{
  
  /**
   * Transforms even the sequence 0,1,2,3,... into reasonably good random numbers
   */
  inline unsigned int get_random_uint(unsigned int seed)
  {
    unsigned int i=(seed^0xA3C59AC3u)*2654435769u;
    i^=(i>>16);
    i*=2654435769u;
    i^=(i>>16);
    i*=2654435769u;
    return i;
  }
  
  
  inline float get_random_float(unsigned int seed)
  {
    return get_random_uint(seed)/(float)UINT_MAX;
  }
  
  
  inline float get_random_float(unsigned int seed, float lower, float upper)
  {
    return ( (upper-lower)*get_random_uint(seed)/(float)UINT_MAX + lower);
  }
  
} // end namespace details


/**
 * Initialize Computational Domain.
 * This function sets up the computational grid and fills
 * liquid particles into the liquid phase.
 *
 * @param grid_width        The length and height size of the grid.
 * @param grid_resolution   The number of grid cells in x and y-directions
 * @param phi               A signed distance field representing the sold boundaries.
 * @param occupied          A signed distance field representing the region to be filled with water.
 */
template<typename F1, typename F2>
inline void initialize_scene(
                             Simulator & simulator
                             , float const & grid_width
                             , int const & grid_resolution
                             , DistanceFieldExp<F1> const & phi
                             , DistanceFieldExp<F2> const & occupied
                             )
{
  initialize_grids(simulator, grid_width, grid_resolution, grid_resolution);
  
  initialize_solid_walls(simulator, phi );
  
  float const wall_threshold = simulator.m_particle_radius*simulator.m_push_fraction;
  
  // Stick some liquid particles in the domain, approx 20 particles inside a fluid cell
  int const K = grid_resolution*grid_resolution*20;
  
  for(int i = 0; i < K; ++i)
  {
    float const x = details::get_random_float(i*2, 0, 1);
    float const y = details::get_random_float(i*2+1, 0, 1);
    
    Vec2f const pt(x,y);
    
    bool const outside_walls = phi(pt) > wall_threshold;
    bool const inside_liquid = occupied(pt) > 0.0;
    
    if(  outside_walls && inside_liquid )
      simulator.m_particles.push_back(pt);
  }
}

// SIMULATOR_INITIALIZE_SCENE_H
#endif
