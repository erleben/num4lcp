#ifndef SIMULATOR_DRAW_H
#define SIMULATOR_DRAW_H

#include <grid2d_draw_zero_levelset.h>
#include <grid2d_draw_grid_lines.h>
#include <grid2d_compute_cell_center_position.h>

#include <simulator.h>
#include <simulator_interpolate_velocity.h>

#include <util_svg_file.h>


inline void draw(Simulator const & sim
                          , bool draw_grid
                          , bool draw_walls
                          , bool draw_liquid
                          , bool draw_particles
                          , bool draw_velocities
                          , bool draw_unilateral
                          , bool draw_solution
                          , bool draw_divergence
                          , SVGFile * svg
                          )
{
  typedef Vec<2,float> V;

  if(draw_grid)
  {
    svg->set_color(0.1, 0.1, 0.5);
    svg->set_stroke_width(1.0);
    draw_grid_lines(sim.m_phi_solid, svg);

    svg->set_color(0.5, 0.1, 0.1);
    svg->set_stroke_width(1.0);
    draw_grid_lines(sim.m_phi_liquid, svg);
  }

  if(draw_walls)
  {
    svg->set_color(0.01, 0.5, 0.01);
    svg->set_stroke_width(2.0);
    draw_zero_levelset( sim.m_phi_solid, svg );
  }

  if(draw_liquid)
  {
    svg->set_color(0.7, 0.1, 0.1);
    draw_zero_levelset( sim.m_phi_liquid, svg );
  }

  if(draw_particles)
  {
    svg->set_color(0.1, 0.1, 0.7);
    for(unsigned int p = 0; p < sim.m_particles.size(); ++p)
    {
      svg->draw_circle(sim.m_particles[p][0] , sim.m_particles[p][1], sim.m_particle_radius);
    }
  }

  if(draw_velocities)
  {
    svg->set_color(0.7, 0.1, 0.7);

    for(int j = 0;j < sim.m_phi_liquid.J(); ++j)
      for(int i = 0; i < sim.m_phi_liquid.I(); ++i)
      {
        V      const start          = compute_cell_center_position(sim.m_phi_liquid,i,j);
        V      const end            = start + 0.01f * interpolate_velocity(sim, start);
        double const arrow_head_len = 0.01f*sim.m_dx;

        svg->draw_arrow2d(start[0], start[1], end[0], end[1], arrow_head_len);

      }
  }

  if(draw_unilateral)
  {
    if(sim.m_is_bilateral.size()>0)
    {
      svg->set_color(0.7, 0.1, 0.1);

      for(int j = 0;j < sim.m_phi_liquid.J(); ++j)
        for(int i = 0; i < sim.m_phi_liquid.I(); ++i)
        {
          V const p = compute_cell_center_position(sim.m_phi_liquid,i,j);

          unsigned int index = i + sim.m_I*j;

          if(!sim.m_is_bilateral[index])
          {
            svg->draw_circle(p[0] , p[1], sim.m_dx/2.0f);
          }
        }
    }
  }


  if(draw_solution)
  {
    if(sim.m_x.size()>0)// Debug visualization???
    {
      float min_p = 1000000.0f;
      float max_p = -min_p;
      for(int j = 0;j < sim.m_phi_liquid.J(); ++j)
        for(int i = 0; i < sim.m_phi_liquid.I(); ++i)
        {
          unsigned int const index = i + sim.m_I*j;
          float const p = sim.m_x[index];
          max_p = (p > max_p) ? p: max_p;
          min_p = (p < min_p) ? p: min_p;
        }

      for(int j = 0;j < sim.m_phi_liquid.J(); ++j)
        for(int i = 0; i < sim.m_phi_liquid.I(); ++i)
        {
          V const pos = compute_cell_center_position(sim.m_phi_liquid,i,j);
          unsigned int const index = i + sim.m_I*j;
          double const p = (sim.m_x[index] - min_p) / (max_p - min_p);
          svg->set_color( p, 0.01, 0.01);
          svg->draw_circle(pos[0] , pos[1], sim.m_dx/2.0f);
        }
    }
  }

  if(draw_divergence)
  {
    if(sim.m_b.size()>0)// Debug visualization???
    {
      float min_b = 1000000.0f;
      float max_b = -min_b;
      for(int j = 0;j < sim.m_phi_liquid.J(); ++j)
        for(int i = 0; i < sim.m_phi_liquid.I(); ++i)
        {
          unsigned int const index = i + sim.m_I*j;
          float const b = sim.m_b[index];
          max_b = (b > max_b)? b: max_b;
          min_b = (b < min_b)? b: min_b;
        }
      for(int j = 0;j < sim.m_phi_liquid.J(); ++j)
        for(int i = 0; i < sim.m_phi_liquid.I(); ++i)
        {
          V const pos = compute_cell_center_position(sim.m_phi_liquid,i,j);
          unsigned int const index = i + sim.m_I*j;
          double const b = (sim.m_b[index] - min_b) / (max_b - min_b);
          svg->set_color( 0.01, 0.01, b);
          svg->draw_circle(pos[0] , pos[1], sim.m_dx/2.0f);
        }
    }
  }
  
}

// SIMULATOR_DRAW_H
#endif
