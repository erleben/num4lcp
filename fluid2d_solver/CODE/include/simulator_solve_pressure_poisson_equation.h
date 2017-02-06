#ifndef SIMULATOR_SOLVE_PRESSURE_POISSON_EQUATION_H
#define SIMULATOR_SOLVE_PRESSURE_POISSON_EQUATION_H

#include <grid2d_interpolate_value.h>
#include <grid2d_compute_cell_center_position.h>

#include <simulator.h>
#include <simulator_compute_cell_has_liquid.h>
#include <simulator_compute_cell_has_empty_space.h>

#include <util_params.h>
#include <util_fraction_inside.h>
#include <util_profiling.h>
#include <util_matlab_write_matrix.h>

#include <solver_profiling_data.h>
#include <solver_solve_eq.h>
#include <solver_solve_lcp.h>

#include <cassert>
#include <iostream>
#include <cstdio>


inline void allocate_system(Simulator & sim)
{
  sim.m_N = sim.m_I*sim.m_J;
  sim.m_K = 0u;
  sim.m_row_indices.resize(sim.m_N*5u,0u);      // At most 5 nonzeros per row in A
  sim.m_column_indices.resize(sim.m_N*5u,0u);
  sim.m_values.resize(sim.m_N*5u,0.0f);

  sim.m_x.resize(sim.m_N,0.0f);
  sim.m_b.resize(sim.m_N,0.0f);
  sim.m_is_bilateral.resize(sim.m_N,false);
}


inline void add_element(Simulator & S, unsigned int i, unsigned int j, float a)
{
  S.m_row_indices[S.m_K]    = i;
  S.m_column_indices[S.m_K] = j;
  S.m_values[S.m_K]         = a;
  S.m_K                     = S.m_K + 1u;
}


inline void system_assembly(Simulator & sim, float const & dt)
{
  typedef typename Simulator::V V;

  using std::max;

  assert( sim.m_K == 0u                || !"system_assembly(): Internal error");
  assert( sim.m_N == (sim.m_I*sim.m_J) || !"system_assembly(): Internal error");

  // Fix bottom boundary
  {
    int const j = 0;
    for(int i = 0; i < sim.m_I; ++i)
    {
      int const index = i + sim.m_I*j;
      add_element(sim, index, index , 1.0f);
      sim.m_is_bilateral[index] = true;
      sim.m_b[index] = 0;
      sim.m_x[index] = 0;
    }
  }
  // Fix left boundary
  {
    int const i = 0;
    for(int j = 1; j < sim.m_J; ++j)
    {
      int const index = i + sim.m_I*j;
      add_element(sim, index, index , 1.0f);
      sim.m_is_bilateral[index] = true;
      sim.m_b[index] = 0;
      sim.m_x[index] = 0;
    }
  }
  // Fix right boundary
  {
    int const i = sim.m_I-1;
    for(int j = 1; j < sim.m_J; ++j)
    {
      int const index = i + sim.m_I*j;
      add_element(sim, index, index , 1.0f);
      sim.m_is_bilateral[index] = true;
      sim.m_b[index] = 0;
      sim.m_x[index] = 0;
    }
  }
  // Fix top boundary
  {
    int const j = sim.m_J-1;
    for(int i = 1; i < sim.m_I-1; ++i)
    {
      int const index = i + sim.m_I*j;
      add_element(sim, index, index , 1.0f);
      sim.m_is_bilateral[index] = true;
      sim.m_b[index] = 0;
      sim.m_x[index] = 0;
    }
  }
  // Fill in interior domain
  for(int j = 1; j < sim.m_J-1; ++j)
  {
    for(int i = 1; i < sim.m_I-1; ++i)
    {
      int const index = i + sim.m_I*j;
      sim.m_b[index] = 0;
      sim.m_x[index] = 0;
      sim.m_is_bilateral[index] = false;

      bool const has_empty_space = compute_cell_has_empty_space(sim, i, j);
      bool const has_liquid      = compute_cell_has_liquid(sim, i, j);

      if( has_empty_space & has_liquid )
      {
        bool const right_has_liquid = compute_cell_has_liquid(sim, i+1,   j);
        bool const left_has_liquid  = compute_cell_has_liquid(sim, i-1,   j);
        bool const up_has_liquid    = compute_cell_has_liquid(sim,   i, j+1);
        bool const down_has_liquid  = compute_cell_has_liquid(sim,   i, j-1);

        float const phi_center = sim.m_phi_liquid(i,j);     // m_phi_liquid is colocated with pressure field
        float const phi_right  = sim.m_phi_liquid(i+1,j);
        float const phi_left   = sim.m_phi_liquid(i-1,j);
        float const phi_up     = sim.m_phi_liquid(i,j+1);
        float const phi_down   = sim.m_phi_liquid(i,j-1);

        float const liquid_fraction_right = fraction_inside(phi_center, phi_right);
        float const liquid_fraction_left  = fraction_inside(phi_center, phi_left);
        float const liquid_fraction_up    = fraction_inside(phi_center, phi_up);
        float const liquid_fraction_down  = fraction_inside(phi_center, phi_down);

        float const b_right = - sim.m_u_weights(i+1,  j) * sim.m_u(i+1,  j) / sim.m_dx;
        float const b_left  =   sim.m_u_weights(  i,  j) * sim.m_u(  i,  j) / sim.m_dx;
        float const b_up    = - sim.m_v_weights(  i,j+1) * sim.m_v(  i,j+1) / sim.m_dx;
        float const b_down  =   sim.m_v_weights(  i,  j) * sim.m_v(  i,  j) / sim.m_dx;

        float const M_right =   sim.m_u_weights(i+1,  j) * dt / sqr(sim.m_dx);
        float const M_left  =   sim.m_u_weights(  i,  j) * dt / sqr(sim.m_dx);
        float const M_up    =   sim.m_v_weights(  i,j+1) * dt / sqr(sim.m_dx);
        float const M_down  =   sim.m_v_weights(  i,  j) * dt / sqr(sim.m_dx);

        float const A_center_right = ((phi_right < 0.0f) ?  M_right : (M_right / max(liquid_fraction_right, 0.01f) ));
        float const A_center_left  = ((phi_left  < 0.0f) ?  M_left  : (M_left  / max(liquid_fraction_left,  0.01f) ));
        float const A_center_up    = ((phi_up    < 0.0f) ?  M_up    : (M_up    / max(liquid_fraction_up,    0.01f) ));
        float const A_center_down  = ((phi_down  < 0.0f) ?  M_down  : (M_down  / max(liquid_fraction_down,  0.01f) ));

        float const A_center       = A_center_right + A_center_left + A_center_up + A_center_down;

        assert( A_center > 0 || !"system_assembly: fatal error");

        float const A_right        = right_has_liquid ? - M_right : 0.0f;
        float const A_left         = left_has_liquid  ? - M_left  : 0.0f;
        float const A_up           = up_has_liquid    ? - M_up    : 0.0f;
        float const A_down         = down_has_liquid  ? - M_down  : 0.0f;

        sim.m_b[index] = b_right + b_left + b_up + b_down;

        add_element(sim,index, index - sim.m_I, A_down  );
        add_element(sim,index, index - 1      , A_left  );
        add_element(sim,index, index          , A_center);
        add_element(sim,index, index + 1      , A_right );
        add_element(sim,index, index + sim.m_I, A_up    );

        // Kenny: Interior fluid cells are forced to be bilateral cells
        // only cells close to a wall or free surface are allowed to be
        // unilateral:
        //
        //    bool const deep_inside = phi_center < -sim.m_dx;
        //
        // However phi_liquid is not a very accurate signed distance field
        // so we need so do something else

        V     const p             = compute_cell_center_position(sim.m_phi_liquid,i,j);
        float const phi_solid     = interpolate_value( sim.m_phi_solid, p );
        bool const  deep_inside    = phi_solid > 1.0*sim.m_dx;
        sim.m_is_bilateral[index] = deep_inside ? true : false;

      }
      else
      {
        add_element(sim,index, index , 1.0f);
        sim.m_is_bilateral[index] = true;
      }

    }
  }
}

inline void eq_solve(Simulator & sim, Params const & params)
{
  if(sim.m_N == 0 || sim.m_K == 0 )
    return;

  ProfilingData profiling_data;

  if(params.host())
  {

    solve_eq_host(
                  sim.m_N
                  , sim.m_K
                  , sim.m_row_indices
                  , sim.m_column_indices
                  , sim.m_values
                  , sim.m_x
                  , sim.m_b
                  , params
                  , profiling_data
                  );

  }else{

    solve_eq_device(
                    sim.m_N
                    , sim.m_K
                    , sim.m_row_indices
                    , sim.m_column_indices
                    , sim.m_values
                    , sim.m_x
                    , sim.m_b
                    , params
                    , profiling_data
                    );

  }

  if(params.profiling() && params.record_time())
  {
    RECORD("host_initialization", profiling_data.m_host_initialization_time );
    RECORD("host_to_device",      profiling_data.m_host_to_device_time      );
    RECORD("device_to_host",      profiling_data.m_device_to_host_time      );
    RECORD("computation",         profiling_data.m_computation_time         );
    RECORD("eq_status",           profiling_data.m_eq_status                );
    RECORD("eq_iteration",        profiling_data.m_eq_iteration             );
  }
  if(params.profiling() && params.record_convergence())
  {
    std::vector<float> tmp = make_std_vector(profiling_data.m_eq_residuals);
    RECORD_VECTOR("eq_residual",  tmp);
  }
}

inline void lcp_solve(Simulator & sim, Params const & params)
{
  ProfilingData profiling_data;

  if(params.host())
  {

    solve_lcp_host( sim.m_N
                   , sim.m_K
                   , sim.m_row_indices
                   , sim.m_column_indices
                   , sim.m_values
                   , sim.m_x
                   , sim.m_b
                   , sim.m_is_bilateral
                   , params
                   , profiling_data
                   );

  }else{
    solve_lcp_device( sim.m_N
                     , sim.m_K
                     , sim.m_row_indices
                     , sim.m_column_indices
                     , sim.m_values
                     , sim.m_x
                     , sim.m_b
                     , sim.m_is_bilateral
                     , params
                     , profiling_data
                     );

  }

  if(params.profiling() && params.record_time())
  {
    RECORD("host_initialization", profiling_data.m_host_initialization_time );
    RECORD("host_to_device",      profiling_data.m_host_to_device_time        );
    RECORD("device_to_host",      profiling_data.m_device_to_host_time        );
    RECORD("computation",         profiling_data.m_computation_time           );
    RECORD("warmstart",           profiling_data.m_psor_time                  );
    RECORD("newton_status",       profiling_data.m_newton_status            );
    RECORD("newton_iteration",    profiling_data.m_newton_iteration         );
    RECORD("psor_status",         profiling_data.m_psor_status              );
    RECORD("psor_iteration",      profiling_data.m_psor_iteration           );

  }
  if(params.profiling() && params.record_convergence())
  {
    std::vector<float> tmp = make_std_vector(profiling_data.m_newton_residuals);
    RECORD_VECTOR("newton_residual", tmp);
    tmp = make_std_vector(profiling_data.m_psor_residuals);
    RECORD_VECTOR("psor_residual", tmp );
  }
}

inline void velocity_update(Simulator & sim, float const & dt)
{
  std::fill(sim.m_u_valid.begin(), sim.m_u_valid.end(), 0 );

  for(int j = 0; j < sim.m_J; ++j)
    for(int i = 1; i < sim.m_I-1; ++i)
    {
      int   const index     = i + j*sim.m_I;

      bool const right_has_liquid = compute_cell_has_liquid(sim, i,   j);
      bool const left_has_liquid  = compute_cell_has_liquid(sim, i-1,   j);

      float const phi_left  = sim.m_phi_liquid(i-1,j);
      float const phi_right = sim.m_phi_liquid(i,  j);

      if(sim.m_u_weights(i,j) > 0 && (right_has_liquid || left_has_liquid))
      {
        bool  const partial_filled = ( phi_left >= 0 || phi_right >= 0);
        float       theta          = partial_filled ? fraction_inside(phi_left, phi_right) : 1.0f;
        theta          = (theta < 0.01f) ? 0.01f : theta;
        float const p_left         = sim.m_x[index-1];
        float const p_right        = sim.m_x[index];

        sim.m_u(i,j) -= dt  * (p_right - p_left) / sim.m_dx / theta;
        sim.m_u_valid(i,j) = 1;
      }
      else
      {
        sim.m_u(i,j) = 0;
      }
    }

  std::fill(sim.m_v_valid.begin(), sim.m_v_valid.end(), 0 );

  for(int j = 1; j < sim.m_J-1; ++j)
    for(int i = 0; i < sim.m_I; ++i)
    {
      int   const index     = i + j*sim.m_I;

      bool const up_has_liquid    = compute_cell_has_liquid(sim,   i, j);
      bool const down_has_liquid  = compute_cell_has_liquid(sim,   i, j-1);

      float const phi_down  = sim.m_phi_liquid(i,j-1);
      float const phi_up    = sim.m_phi_liquid(i,  j);

      if(sim.m_v_weights(i,j) > 0 && (up_has_liquid || down_has_liquid))
      {
        bool  const partial_filled = ( phi_down >= 0 || phi_up >= 0);
        float       theta          = partial_filled ? fraction_inside(phi_down, phi_up) : 1.0f;
        theta          = (theta < 0.01f) ? 0.01f : theta;
        float const p_down         = sim.m_x[index-sim.m_I];
        float const p_up           = sim.m_x[index];

        sim.m_v(i,j) -= dt  * (p_up - p_down) / sim.m_dx / theta;
        sim.m_v_valid(i,j) = 1;
      }
      else
      {
        sim.m_v(i,j) = 0;
      }
    }
}

inline void solve_pressure_poisson_equation(Simulator & sim, float const & dt, Params const & params)
{
  START_TIMER("assembly_allocation");
  allocate_system(sim);
  STOP_TIMER("assembly_allocation");

  START_TIMER("assembly");
  system_assembly(sim,dt);
  STOP_TIMER("assembly");

  std::ofstream file;

  if(params.profiling() && params.record_data())
  {
    static unsigned int file_count = 0u;
    std::stringstream filename;
    filename << params.output_path() << params.prefix() << "data" << params.postfix() << "_";
    filename.fill('0');
    filename.width(4);
    filename << ++file_count;
    filename << ".m";

    file.open( filename.str().c_str() );

    file << "clear all;" << std::endl;
    file << "close all;" << std::endl;
    file << "clc;"       << std::endl;

    file << "N = " <<  sim.m_I << std::endl;
    file << matlab_write_matrix("A", sim.m_row_indices, sim.m_column_indices, sim.m_values, sim.m_K, sim.m_N, sim.m_N) << std::endl;
    file << "D = A - A';" << std::endl;
    file << "min(D(:))" << std::endl;
    file << "figure(1); spy(A); title('Coefficient Matrix')" << std::endl;

    file << "lambda = sort( eig( full(A) ) );" << std::endl;
    file << "min_eigen_value = min( lambda(:) )" << std::endl;
    file << "figure(2); plot( lambda ); colormap gray; title('Sorted Eigenvalues')" << std::endl;


    file << matlab_write_vector("b", sim.m_b ) << std::endl;
    file << "B = reshape(b,N,N)';" << std::endl;
    file << "figure(3); imagesc(B); colorbar; colormap gray; title('Divergence field')" << std::endl;

    file << matlab_write_vector("s", sim.m_is_bilateral ) << std::endl;
    file << "S = reshape(s,N,N)';" << std::endl;
    file << "figure(4); imagesc(S); colorbar; colormap gray; title('Bilateral constraints')" << std::endl;

  }

  START_TIMER("pressure_solver");
  if(params.use_lcp())
  {
    lcp_solve(sim, params);
  }
  else
  {
    eq_solve(sim, params);
  }
  STOP_TIMER("pressure_solver");

  if(params.profiling() && params.record_data())
  {
    file << matlab_write_vector("x", sim.m_x ) << std::endl;
    file << "X = reshape(x,N,N)';" << std::endl;
    file << "figure(5); imagesc(X); colorbar; colormap gray; title('Iterative solution')" << std::endl;

    file << "p = A / b;" << std::endl;
    file << "P = reshape(p,N,N)';" << std::endl;
    file << "figure(6); imagesc(P); colorbar; colormap gray; title('Accurate solution')" << std::endl;
    
    file.flush();
    file.close();
  }
  
  START_TIMER("velocity_update");
  velocity_update(sim,dt);
  STOP_TIMER("velocity_update");
}

// SIMULATOR_SOLVE_PRESSURE_POISSON_EQUATION_H
#endif
