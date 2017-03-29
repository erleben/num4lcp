#ifndef UTIL_PARAMS_H
#define UTIL_PARAMS_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>  // needed for std::max

namespace detail
{
  inline void get_strings_from_text_list(std::string const & text_list, std::vector<std::string> & tokens)
  {
    std::stringstream ss(text_list);

    std::string token;

    tokens.clear();

    while (std::getline(ss,token, ' '))
    {
      tokens.push_back(token);
      std::cout << "token = " << token << std::endl;
    }
  }

  inline bool string2bool (std::string const & value)
  {
    return !value.empty() &&
    (value.compare("true") == 0 ||
     atoi (value.c_str ()) != 0);
    // ^std=c++0x: std::stoi(value) != 0
  }

  inline std::string bool2string (bool const & value)
  {
    return (value ? "true" : "false");
  }

  template<typename VAL>
  inline VAL string2number(std::string const & value)
  {
    VAL result;

    std::stringstream stream(value);
    stream >> result;
    if (stream.fail())
    {
      std::cerr << "string2number(): could not convert string to a number" << std::endl;
    }
    return result;
  }
}// namespace detail


class Params
{
private:

  typedef float        T;
  typedef unsigned int I;

protected:

  T             m_total_time;                     ///< The total simulation time after which the application stop simulating and closes down. If zero
                                                  ///< the application runs until the end-user exits the application. Negative values are clamped to zero.
  T             m_time_step;                      ///< The time step size  used in case fixed time step size is enabled. Negative values
                                                  ///< are clamped to default time-step value of 0.01 seconds
  
  
  bool          m_host;                            ///< Boolean flag indicating whether we are running on host (CPU) or device (GPU).
  bool          m_use_lcp;                         ///< Boolean flag that flips whether slip BCs or separating BCs are used (linear system versus linear complementarity problem).
  I             m_device_number;                   ///< Unsigned integer used to select which GPU device the computation should be performed.

  I             m_preconditioner_choice;           ///< Unsigned integer used to select what type of preconditioner to use.
  
  I             m_newton_max_iterations;           ///< The maximum allowed outer solver iterations. This is the maximum allowed number of Newton iterations.
  T             m_newton_absolute_tolerance;       ///< The absolute tolerance level for the Newton solver
  T             m_newton_relative_tolerance;       ///< The relative tolerance level for the Newton solver
  T             m_newton_stagnation_tolerance;     ///< The stagnation tolerance level for the Newton solver
  T             m_newton_rho;                      ///< The relative residual stop criteria for sub solver used inside the Newton solver.

  I             m_sub_solver_choice;               ///< Unsigned integer used to select what type of subsolsver to use.
  I             m_sub_solver_max_iterations;       ///< The maximum allowed sub solver iterations. This is the maximum allowed number of PCG iterations.
  T             m_sub_solver_absolute_tolerance;   ///< The absolute tolerance level for the PCG solver.
  T             m_sub_solver_relative_tolerance;   ///< The relative tolerance level for the PCG solver.

  I             m_line_search_max_iterations;      ///< Maximum number of itreations allowed in line search method
  T             m_line_search_alpha;               ///< The step reduction parameter
  T             m_line_search_beta;                ///< The sufficient descrease coefficient
  T             m_line_search_gamma;               ///< To small squared step length guard threshold value.
  
  bool          m_psor_warmstart;                  ///< Boolean flag that indicates if pgs warmstart should be used for LCP solving
  I             m_psor_max_iterations;             ///< The maximum allowed PSOR iterations.
  T             m_psor_absolute_tolerance;         ///< The absolute tolerance level for the PSOR solver
  T             m_psor_relative_tolerance;         ///< The relative tolerance level for the PSOR solver
  T             m_psor_stagnation_tolerance;       ///< The stagnation tolerance level for the PSOR solver
  T             m_psor_relaxation;                 ///< The PSOR relaxation.

  bool          m_grad_warmstart;                  ///< Boolean flag that indicates if sd warmstart should be used for LCP solving
  I             m_grad_max_iterations;             ///< The maximum allowed sd iterations.
  T             m_grad_absolute_tolerance;         ///< The absolute tolerance level for the sd solver
  T             m_grad_relative_tolerance;         ///< The relative tolerance level for the sd solver
  T             m_grad_stagnation_tolerance;       ///< The stagnation tolerance level for the sd solver

  I             m_eq_max_iterations;               ///< The maximum allowed sub solver iterations. This is the maximum allowed number of iterations.
  T             m_eq_absolute_tolerance;           ///< The absolute tolerance level for the solver.
  T             m_eq_relative_tolerance;           ///< The relative tolerance level for the solver.
  
  
  std::string   m_scene;                           ///< The string indicating which scene to create
  I             m_grid_resolution;                 ///< The grid resolution of the computation domain.
  T             m_grid_length;                     ///< The length fo the grid side.
  
  bool          m_profiling;                       ///< Enables profiling of simulator
  bool          m_record_data;                     ///< When profiling is enabled this enables recording pressure solve data
  bool          m_record_time;                     ///< When profiling is enabled this enables recording timing results
  bool          m_record_convergence;              ///< When profiling is enabled this enables recording convergence results
  std::string   m_output_path;                     ///< Output path for profiling data and other data files etc.
  
  bool          m_draw_grid;
  bool          m_draw_walls;
  bool          m_draw_liquid;
  bool          m_draw_particles;
  bool          m_draw_velocities;
  bool          m_draw_unilateral;
  bool          m_draw_solution;
  bool          m_draw_divergence;

  bool          m_verbose;                        ///< Boolean flag used to toggle output messages on and off.

  std::string   m_prefix;
  std::string   m_postfix;
  
public:

  T           const & total_time()                    const { return this->m_total_time;                    }
  T           const & time_step()                     const { return this->m_time_step;                     }
  bool        const & host()                          const { return this->m_host;                          }
  bool                use_gpu()                       const { return !(this->m_host);                       }
  bool        const & use_lcp()                       const { return this->m_use_lcp;                       }
  I           const & device_number()                 const { return this->m_device_number;                 }
  
  I           const & preconditioner_choice()         const { return this->m_preconditioner_choice;         }

  I           const & newton_max_iterations()         const { return this->m_newton_max_iterations;         }
  T           const & newton_absolute_tolerance()     const { return this->m_newton_absolute_tolerance;     }
  T           const & newton_relative_tolerance()     const { return this->m_newton_relative_tolerance;     }
  T           const & newton_stagnation_tolerance()   const { return this->m_newton_stagnation_tolerance;   }
  T           const & newton_rho()                    const { return this->m_newton_rho;                    }
  
  I           const & sub_solver_choice()             const { return this->m_sub_solver_choice;             }
  I           const & sub_solver_max_iterations()     const { return this->m_sub_solver_max_iterations;     }
  T           const & sub_solver_absolute_tolerance() const { return this->m_sub_solver_absolute_tolerance; }
  T           const & sub_solver_relative_tolerance() const { return this->m_sub_solver_relative_tolerance; }

  I           const & line_search_max_iterations()    const { return this->m_line_search_max_iterations;    }
  T           const & line_search_alpha()             const { return this->m_line_search_alpha;             }
  T           const & line_search_beta()              const { return this->m_line_search_beta;              }
  T           const & line_search_gamma()             const { return this->m_line_search_gamma;             }

  bool        const & psor_warmstart()                const { return this->m_psor_warmstart;                }
  I           const & psor_max_iterations()           const { return this->m_psor_max_iterations;           }
  T           const & psor_absolute_tolerance()       const { return this->m_psor_absolute_tolerance;       }
  T           const & psor_relative_tolerance()       const { return this->m_psor_relative_tolerance;       }
  T           const & psor_stagnation_tolerance()     const { return this->m_psor_stagnation_tolerance;     }
  T           const & psor_relaxation()               const { return this->m_psor_relaxation;               }

  bool        const & grad_warmstart()                const { return this->m_grad_warmstart;                }
  I           const & grad_max_iterations()           const { return this->m_grad_max_iterations;           }
  T           const & grad_absolute_tolerance()       const { return this->m_grad_absolute_tolerance;       }
  T           const & grad_relative_tolerance()       const { return this->m_grad_relative_tolerance;       }
  T           const & grad_stagnation_tolerance()     const { return this->m_grad_stagnation_tolerance;     }

  I           const & eq_max_iterations()             const { return this->m_eq_max_iterations;             }
  T           const & eq_absolute_tolerance()         const { return this->m_eq_absolute_tolerance;         }
  T           const & eq_relative_tolerance()         const { return this->m_eq_relative_tolerance;         }
  
  std::string const & scene()                         const { return this->m_scene;                         }
  I           const & grid_resolution()               const { return this->m_grid_resolution;               }
  T           const & grid_length()                   const { return this->m_grid_length;                   }
  
  bool        const & profiling()                     const { return this->m_profiling;                     }
  bool        const & record_data()                   const { return this->m_record_data;                   }
  bool        const & record_time()                   const { return this->m_record_time;                   }
  bool        const & record_convergence()            const { return this->m_record_convergence;            }
  std::string const & output_path()                   const { return this->m_output_path;                   }
  
  bool        const & draw_grid()                     const { return this->m_draw_grid;                     }
  bool        const & draw_walls()                    const { return this->m_draw_walls;                    }
  bool        const & draw_liquid()                   const { return this->m_draw_liquid;                   }
  bool        const & draw_particles()                const { return this->m_draw_particles;                }
  bool        const & draw_velocities()               const { return this->m_draw_velocities;               }
  bool        const & draw_unilateral()               const { return this->m_draw_unilateral;               }
  bool        const & draw_solution()                 const { return this->m_draw_solution;                 }
  bool        const & draw_divergence()               const { return this->m_draw_divergence;               }
  bool        const & verbose()                       const { return this->m_verbose;                       }

  std::string const & prefix()                        const { return this->m_prefix;                        }
  std::string const & postfix()                       const { return this->m_postfix;                       }

public:
  
  T            & total_time()                     { return this->m_total_time;                    }
  T            & time_step()                      { return this->m_time_step;                     }
  bool         & host()                           { return this->m_host;                          }
  bool         & use_lcp()                        { return this->m_use_lcp;                       }
  I            & device_number()                  { return this->m_device_number;                 }
  
  I            & preconditioner_choice()          { return this->m_preconditioner_choice;         }
  
  I            & newton_max_iterations()          { return this->m_newton_max_iterations;         }
  T            & newton_absolute_tolerance()      { return this->m_newton_absolute_tolerance;     }
  T            & newton_relative_tolerance()      { return this->m_newton_relative_tolerance;     }
  T            & newton_stagnation_tolerance()    { return this->m_newton_stagnation_tolerance;   }
  T            & newton_rho()                     { return this->m_newton_rho;                    }
  
  I            & sub_solver_choice()              { return this->m_sub_solver_choice;             }
  I            & sub_solver_max_iterations()      { return this->m_sub_solver_max_iterations;     }
  T            & sub_solver_absolute_tolerance()  { return this->m_sub_solver_absolute_tolerance; }
  T            & sub_solver_relative_tolerance()  { return this->m_sub_solver_relative_tolerance; }
  
  I            & line_search_max_iterations()     { return this->m_line_search_max_iterations;    }
  T            & line_search_alpha()              { return this->m_line_search_alpha;             }
  T            & line_search_beta()               { return this->m_line_search_beta;              }
  T            & line_search_gamma()              { return this->m_line_search_gamma;             }
  
  bool         & psor_warmstart()                 { return this->m_psor_warmstart;                }
  I            & psor_max_iterations()            { return this->m_psor_max_iterations;           }
  T            & psor_absolute_tolerance()        { return this->m_psor_absolute_tolerance;       }
  T            & psor_relative_tolerance()        { return this->m_psor_relative_tolerance;       }
  T            & psor_stagnation_tolerance()      { return this->m_psor_stagnation_tolerance;     }
  T            & psor_relaxation()                { return this->m_psor_relaxation;               }

  bool         & grad_warmstart()                 { return this->m_grad_warmstart;                }
  I            & grad_max_iterations()            { return this->m_grad_max_iterations;           }
  T            & grad_absolute_tolerance()        { return this->m_grad_absolute_tolerance;       }
  T            & grad_relative_tolerance()        { return this->m_grad_relative_tolerance;       }
  T            & grad_stagnation_tolerance()      { return this->m_grad_stagnation_tolerance;     }


  I            & eq_max_iterations()              { return this->m_eq_max_iterations;             }
  T            & eq_absolute_tolerance()          { return this->m_eq_absolute_tolerance;         }
  T            & eq_relative_tolerance()          { return this->m_eq_relative_tolerance;         }
  
  std::string  & scene()                          { return this->m_scene;                         }
  I            & grid_resolution()                { return this->m_grid_resolution;               }
  T            & grid_length()                    { return this->m_grid_length;                   }
  
  bool         & profiling()                      { return this->m_profiling;                     }
  bool         & record_data()                    { return this->m_record_data;                   }
  bool         & record_time()                    { return this->m_record_time;                   }
  bool         & record_convergence()             { return this->m_record_convergence;            }
  std::string  & output_path()                    { return this->m_output_path;                   }
  
  bool         & draw_grid()                      { return this->m_draw_grid;                     }
  bool         & draw_walls()                     { return this->m_draw_walls;                    }
  bool         & draw_liquid()                    { return this->m_draw_liquid;                   }
  bool         & draw_particles()                 { return this->m_draw_particles;                }
  bool         & draw_velocities()                { return this->m_draw_velocities;               }
  bool         & draw_unilateral()                { return this->m_draw_unilateral;               }
  bool         & draw_solution()                  { return this->m_draw_solution;                 }
  bool         & draw_divergence()                { return this->m_draw_divergence;               }
  bool         & verbose()                        { return this->m_verbose;                       }

public:

  Params()
  : m_total_time(10.0)
  , m_time_step(0.01)
  , m_host(true)
  , m_use_lcp(false)
  , m_device_number(0u)
  , m_preconditioner_choice(0u)
  , m_newton_max_iterations(10u)
  , m_newton_absolute_tolerance(1e-6)
  , m_newton_relative_tolerance(1e-3)
  , m_newton_stagnation_tolerance(1e-6)
  , m_newton_rho(0.01)
  , m_sub_solver_choice(0u)
  , m_sub_solver_max_iterations(100u)
  , m_sub_solver_absolute_tolerance(1e-6)
  , m_sub_solver_relative_tolerance(1e-3)
  , m_line_search_max_iterations(100u)
  , m_line_search_alpha(0.5)
  , m_line_search_beta(0.001)
  , m_line_search_gamma(0.000001)
  , m_psor_warmstart(false)
  , m_psor_max_iterations(10u)
  , m_psor_absolute_tolerance(1e-6)
  , m_psor_relative_tolerance(1e-3)
  , m_psor_stagnation_tolerance(1e-6)
  , m_psor_relaxation(0.75)
  , m_grad_warmstart(false)
  , m_grad_max_iterations(10u)
  , m_grad_absolute_tolerance(1e-6)
  , m_grad_relative_tolerance(1e-4)
  , m_grad_stagnation_tolerance(1e-6)
  , m_eq_max_iterations(100u)
  , m_eq_absolute_tolerance(1e-6)
  , m_eq_relative_tolerance(1e-3)
  , m_scene("1")
  , m_grid_resolution(32)
  , m_grid_length(1.0)
  , m_profiling(false)
  , m_record_data(false)
  , m_record_time(false)
  , m_record_convergence(false)
  , m_output_path("./")
  , m_draw_grid(true)
  , m_draw_walls(true)
  , m_draw_liquid(true)
  , m_draw_particles(true)
  , m_draw_velocities(true)
  , m_draw_unilateral(false)
  , m_draw_solution(false)
  , m_draw_divergence(false)
  , m_verbose(true)
  , m_prefix("")
  , m_postfix("")
  {
  }
  
  virtual ~Params(){}

private:

  std::string get_string_value(const char * file_name, std::string const & token, std::string const & default_value = "")
  {
    std::ifstream file_stream;
    file_stream.open(file_name);
    std::string whitespaces (" \t\f\v\n\r");

    while (!file_stream.eof())
    {
      std::string line;
      std::getline(file_stream, line);
      size_t delimiter_pos = line.find_first_of(whitespaces);
      std::string token_str = line.substr(0,delimiter_pos);
      if (token_str.compare(token) == 0)
      {
        std::string sub_str = line.substr(line.find_first_of("=")+1);
        std::string param = sub_str.substr(sub_str.find_first_not_of(whitespaces));
        return param.erase(param.find_last_not_of(whitespaces)+1);
      }
    }
    
    return default_value;
  }

public:

  void load(const char *  file_name)
  {
    using std::min;
    using std::max;

    std::string config_file = file_name;

    if( config_file.empty() )
    {
      std::cout << "Params::load() No config file is given - I am using default.cfg" << std::endl;

      config_file   = "default.cfg";
      file_name     = "default.cfg";
    }

    std::ifstream file_stream;
    file_stream.open(file_name);

    if (!file_stream.good())
    {
      std::cout <<  "Params::load(): file not found: " << file_name <<  std::endl;
      return;
    }
    
    std::string total_time                    = get_string_value(file_name, "total_time", "10.0"                    );
    std::string time_step                     = get_string_value(file_name, "time_step", "0.01"                     );
    m_total_time                              = detail::string2number<T>( total_time );
    m_time_step                               = detail::string2number<T>( time_step  );
    
    std::string host                          = get_string_value(file_name, "host", "true"                          );
    std::string use_lcp                       = get_string_value(file_name, "use_lcp", "false"                      );
    std::string preconditioner_choice         = get_string_value(file_name, "preconditioner_choice", "0"            );
    std::string device_number                 = get_string_value(file_name, "device_number", "0"                    );
    m_host                                    = detail::string2bool( host );
    m_use_lcp                                 = detail::string2bool( use_lcp );
    m_preconditioner_choice                   = detail::string2number<I>( preconditioner_choice );
    m_device_number                           = detail::string2number<I>( device_number );
    
    std::string newton_max_iterations         = get_string_value(file_name, "newton_max_iterations", "100"          );
    std::string newton_absolute_tolerance     = get_string_value(file_name, "newton_absolute_tolerance", "1e-6"     );
    std::string newton_relative_tolerance     = get_string_value(file_name, "newton_relative_tolerance", "1e-3"     );
    std::string newton_stagnation_tolerance   = get_string_value(file_name, "newton_stagnation_tolerance", "1e-6"   );
    std::string newton_rho                    = get_string_value(file_name, "newton_rho", "0.1"                     );
    m_newton_max_iterations                   = detail::string2number<I>( newton_max_iterations );
    m_newton_absolute_tolerance               = detail::string2number<T>( newton_absolute_tolerance );
    m_newton_relative_tolerance               = detail::string2number<T>( newton_relative_tolerance );
    m_newton_stagnation_tolerance             = detail::string2number<T>( newton_stagnation_tolerance );
    m_newton_rho                              = detail::string2number<T>( newton_rho );

    std::string sub_solver_choice             = get_string_value(file_name, "sub_solver_choice", "0"                );
    std::string sub_solver_max_iterations     = get_string_value(file_name, "sub_solver_max_iterations", "100"      );
    std::string sub_solver_absolute_tolerance = get_string_value(file_name, "sub_solver_absolute_tolerance", "1e-6" );
    std::string sub_solver_relative_tolerance = get_string_value(file_name, "sub_solver_relative_tolerance", "1e-3" );
    m_sub_solver_choice                       = detail::string2number<I>( sub_solver_choice );
    m_sub_solver_max_iterations               = detail::string2number<I>( sub_solver_max_iterations );
    m_sub_solver_absolute_tolerance           = detail::string2number<T>( sub_solver_absolute_tolerance );
    m_sub_solver_relative_tolerance           = detail::string2number<T>( sub_solver_relative_tolerance );
    
    std::string line_search_max_iterations    = get_string_value(file_name, "line_search_max_iterations", "100"     );
    std::string line_search_alpha             = get_string_value(file_name, "line_search_alpha", "0.5"              );
    std::string line_search_beta              = get_string_value(file_name, "line_search_beta", "1e-4"              );
    std::string line_search_gamma             = get_string_value(file_name, "line_search_gamma", "1e-6"             );
    std::string psor_warmstart                = get_string_value(file_name, "psor_warmstart", "false"               );
    std::string psor_max_iterations           = get_string_value(file_name, "psor_max_iterations", "10"             );
    std::string psor_absolute_tolerance       = get_string_value(file_name, "psor_absolute_tolerance", "1e-6"       );
    std::string psor_relative_tolerance       = get_string_value(file_name, "psor_relative_tolerance", "1e-3"       );
    std::string psor_stagnation_tolerance     = get_string_value(file_name, "psor_stagnation_tolerance", "1e-6"     );
    std::string psor_relaxation               = get_string_value(file_name, "psor_relaxation", "0.75"               );
    std::string grad_warmstart                = get_string_value(file_name, "grad_warmstart", "false"               );
    std::string grad_max_iterations           = get_string_value(file_name, "grad_max_iterations", "10"             );
    std::string grad_absolute_tolerance       = get_string_value(file_name, "grad_absolute_tolerance", "1e-6"       );
    std::string grad_relative_tolerance       = get_string_value(file_name, "grad_relative_tolerance", "1e-3"       );
    std::string grad_stagnation_tolerance     = get_string_value(file_name, "grad_stagnation_tolerance", "1e-6"     );

    m_psor_warmstart                          = detail::string2bool( psor_warmstart );
    m_psor_max_iterations                     = detail::string2number<I>( psor_max_iterations   );
    m_psor_absolute_tolerance                 = detail::string2number<T>( psor_absolute_tolerance    );
    m_psor_relative_tolerance                 = detail::string2number<T>( psor_relative_tolerance    );
    m_psor_stagnation_tolerance               = detail::string2number<T>( psor_stagnation_tolerance   );
    m_psor_relaxation                         = detail::string2number<T>( psor_relaxation );
    m_grad_warmstart                          = detail::string2bool( grad_warmstart );
    m_grad_max_iterations                     = detail::string2number<I>( grad_max_iterations   );
    m_grad_absolute_tolerance                 = detail::string2number<T>( grad_absolute_tolerance    );
    m_grad_relative_tolerance                 = detail::string2number<T>( grad_relative_tolerance    );
    m_grad_stagnation_tolerance               = detail::string2number<T>( grad_stagnation_tolerance   );

    std::string eq_max_iterations             = get_string_value(file_name, "eq_max_iterations", "100"              );
    std::string eq_absolute_tolerance         = get_string_value(file_name, "eq_absolute_tolerance", "1e-6"         );
    std::string eq_relative_tolerance         = get_string_value(file_name, "eq_relative_tolerance", "1e-3"         );
    m_eq_max_iterations                       = detail::string2number<I>( eq_max_iterations );
    m_eq_absolute_tolerance                   = detail::string2number<T>( eq_absolute_tolerance );
    m_eq_relative_tolerance                   = detail::string2number<T>( eq_relative_tolerance );
    
    std::string scene                         = get_string_value(file_name, "scene", "1"                            );
    std::string grid_resolution               = get_string_value(file_name, "grid_resolution", "64"                 );
    std::string grid_length                   = get_string_value(file_name, "grid_length", "1.0"                    );
    m_scene                                   = scene;
    m_grid_resolution                         = detail::string2number<I>( grid_resolution  );
    m_grid_length                             = detail::string2number<T>( grid_length  );
    
    std::string profiling                     = get_string_value(file_name, "profiling", "false"                    );
    std::string record_data                   = get_string_value(file_name, "record_data", "false"                  );
    std::string record_time                   = get_string_value(file_name, "record_time", "false"                  );
    std::string record_convergence            = get_string_value(file_name, "record_convergence", "false"           );
    std::string output_path                   = get_string_value(file_name, "output_path", "./"                     );
    m_profiling                               = detail::string2bool( profiling );
    m_record_data                             = detail::string2bool( record_data );
    m_record_time                             = detail::string2bool( record_time );
    m_record_convergence                      = detail::string2bool( record_convergence );
    m_output_path                             = output_path;

    std::string draw_grid                     = get_string_value(file_name, "draw_grid", "true"                     );
    std::string draw_walls                    = get_string_value(file_name, "draw_walls", "true"                    );
    std::string draw_liquid                   = get_string_value(file_name, "draw_liquid", "true"                   );
    std::string draw_particles                = get_string_value(file_name, "draw_particles", "true"                );
    std::string draw_velocities               = get_string_value(file_name, "draw_velocities", "true"               );
    std::string draw_unilateral               = get_string_value(file_name, "draw_unilateral", "false"              );
    std::string draw_solution                 = get_string_value(file_name, "draw_solution", "false"                );
    std::string draw_divergence               = get_string_value(file_name, "draw_divergence", "false"              );
    std::string verbose                       = get_string_value(file_name, "verbose", "true"                       );
    m_draw_grid                               = detail::string2bool( draw_grid       );
    m_draw_particles                          = detail::string2bool( draw_particles  );
    m_draw_walls                              = detail::string2bool( draw_walls      );
    m_draw_liquid                             = detail::string2bool( draw_liquid     );
    m_draw_velocities                         = detail::string2bool( draw_velocities );
    m_draw_unilateral                         = detail::string2bool( draw_unilateral );
    m_draw_solution                           = detail::string2bool( draw_solution   );
    m_draw_divergence                         = detail::string2bool( draw_divergence );
    m_verbose                                 = detail::string2bool( verbose         );

    m_prefix                                  = get_string_value(file_name, "prefix", ""   );
    m_postfix                                 = get_string_value(file_name, "postfix", ""  );

   }
  
};


inline std::string print_params(Params const & params)
{
  std::stringstream output;

  output << "total_time                  = " << params.total_time()                     << std::endl;
  output << "time_step                   = " << params.time_step()                      << std::endl;
  output << "host                        = " << params.host()                           << std::endl;
  output << "use_lcp                     = " << params.use_lcp()                        << std::endl;
  output << "device_number               = " << params.device_number()                  << std::endl;
  output << "preconditioner_choice       = " << params.preconditioner_choice()          << std::endl;
  output << "# 0: identity, 1: diagonal, 2 : standard scaled bridson, 3: static scaled bridson, 4: scaled bridson Lin strategy" << std::endl;
  output << "newton_max_iterations       = " << params.newton_max_iterations()          << std::endl;
  output << "newton_absolute_tolerance   = " << params.newton_absolute_tolerance()      << std::endl;
  output << "newton_relative_tolerance   = " << params.newton_relative_tolerance()      << std::endl;
  output << "newton_stagnation_tolerance = " << params.newton_stagnation_tolerance()    << std::endl;
  output << "newton_rho                  = " << params.newton_rho()                     << std::endl;
  output << "sub_solver_choice           = " << params.sub_solver_choice()              << std::endl;
  output << "# 0: PCG, 1: CR, 2: BiCG (not supported), 3: BiCGStab, 4: GMRES"           << std::endl;
  output << "sub_solver_max_iterations     = " << params.sub_solver_max_iterations()     << std::endl;
  output << "sub_solver_absolute_tolerance = " << params.sub_solver_absolute_tolerance() << std::endl;
  output << "sub_solver_relative_tolerance = " << params.sub_solver_relative_tolerance() << std::endl;
  output << "line_search_max_iterations    = " << params.line_search_max_iterations()   << std::endl;
  output << "line_search_alpha             = " << params.line_search_alpha()            << std::endl;
  output << "line_search_beta              = " << params.line_search_beta()             << std::endl;
  output << "line_search_gamma             = " << params.line_search_gamma()            << std::endl;
  output << "psor_warmstart                = " << params.psor_warmstart()               << std::endl;
  output << "psor_max_iterations           = " << params.psor_max_iterations()          << std::endl;
  output << "psor_absolute_tolerance       = " << params.psor_absolute_tolerance()      << std::endl;
  output << "psor_relative_tolerenace      = " << params.psor_relative_tolerance()      << std::endl;
  output << "psor_stagnation_tolerance     = " << params.psor_stagnation_tolerance()    << std::endl;
  output << "psor_relaxation               = " << params.psor_relaxation()              << std::endl;
  output << "grad_warmstart                = " << params.grad_warmstart()               << std::endl;
  output << "grad_max_iterations           = " << params.grad_max_iterations()          << std::endl;
  output << "grad_absolute_tolerance       = " << params.grad_absolute_tolerance()      << std::endl;
  output << "grad_relative_tolerenace      = " << params.grad_relative_tolerance()      << std::endl;
  output << "grad_stagnation_tolerance     = " << params.grad_stagnation_tolerance()    << std::endl;
  output << "eq_max_iterations             = " << params.eq_max_iterations()            << std::endl;
  output << "eq_absolute_tolerance         = " << params.eq_absolute_tolerance()        << std::endl;
  output << "eq_relative_tolerance         = " << params.eq_relative_tolerance()        << std::endl;
  output << "scene                         = " << params.scene()                        << std::endl;
  output << "grid_resolution               = " << params.grid_resolution()              << std::endl;
  output << "grid_length                   = " << params.grid_length()                  << std::endl;
  output << "profiling                     = " << params.profiling()                    << std::endl;
  output << "record_data                   = " << params.record_data()                  << std::endl;
  output << "record_time                   = " << params.record_time()                  << std::endl;
  output << "record_convergence            = " << params.record_convergence()           << std::endl;
  output << "output_path                   = " << params.output_path()                  << std::endl;
  output << "draw_grid                     = " << params.draw_grid()                    << std::endl;
  output << "draw_walls                    = " << params.draw_walls()                   << std::endl;
  output << "draw_liquid                   = " << params.draw_liquid()                  << std::endl;
  output << "draw_particles                = " << params.draw_particles()               << std::endl;
  output << "draw_velocities               = " << params.draw_velocities()              << std::endl;
  output << "draw_unilateral               = " << params.draw_unilateral()              << std::endl;
  output << "draw_solution                 = " << params.draw_solution()                << std::endl;
  output << "draw_divergence               = " << params.draw_divergence()              << std::endl;
  output << "verbose                       = " << params.verbose()                      << std::endl;
  output << "prefix                        = " << params.prefix()                       << std::endl;
  output << "postfix                       = " << params.postfix()                      << std::endl;

  output.flush();

  return output.str();

}

// UTIL_PARAMS_H
#endif
