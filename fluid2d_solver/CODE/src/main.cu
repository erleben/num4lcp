#include <gpu_print_device.h>
#include <gpu_print_versions.h>

#include <util_matlab_write_vector.h>
#include <util_matlab_write_profiling.h>
#include <util_distance_field.h>
#include <util_params.h>
#include <util_profiling.h>

#include <simulator.h>
#include <simulator_draw.h>
#include <simulator_initialize_scene.h>
#include <simulator_run.h>

#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

class Application
{
protected:

  Params    m_params;
  Simulator m_simulator;
  float     m_cur_time;

protected:

  void scene1()
  {
    Box     A(0.2f,0.8f,0.4f,0.6f);
    Circle  B(0.5f,0.5f,0.4f);
    Plane   C;
    initialize_scene(m_simulator, m_params.grid_length(), m_params.grid_resolution(),-A%B, C);
  }

  void scene2()
  {
    Circle   A(0.3f,0.5f,0.25f);
    Circle   B(0.7f,0.5f,0.25f);
    initialize_scene(m_simulator, m_params.grid_length(), m_params.grid_resolution(), A+B, -dilation(A+B,-0.0575f));
  }

  void scene3()
  {
    Box     A(0.1f,0.9f,0.1f,0.9f);
    Box     B(0.3f,0.9f,0.6f,0.7f);
    Box     C(0.1f,0.7f,0.4f,0.5f);
    Box     D(0.3f,0.9f,0.2f,0.3f);
    Circle  E(0.8f,0.8f,0.1f);
    initialize_scene(m_simulator, m_params.grid_length(), m_params.grid_resolution(), A % -B % -C % -D , E );
  }

  void scene4()
  {
    Box     A(0.1f,0.9f,0.1f,0.9f);
    Box     B(0.1f,0.3f,0.1f,0.7f);
    initialize_scene(m_simulator, m_params.grid_length(), m_params.grid_resolution(), A , B);
  }

  void scene5()
  {
    Box     A(0.1f,0.9f,0.1f,0.9f);
    Box     B(0.1f,0.3f,0.1f,0.7f);
    Box     C(0.7f,0.9f,0.1f,0.7f);
    initialize_scene(m_simulator, m_params.grid_length(), m_params.grid_resolution(), A , B+C );
  }

public:

  void run(int argc, char **argv)
  {
    if (argv[1])
    {
      std::cout << "Application::run(): Loading cfg-file = " << argv[1] << std::endl;
      m_params.load( argv[1] );
    }
    else
    {
      std::cout << "Application::run(): Loading cfg-file default.cfg " << std::endl;
      m_params.load( "" );
    }

    //--- Setup simulation stuff ----------------
    if (m_params.scene().compare("1")==0 )
    {
      scene1();
    }
    else if (m_params.scene().compare("2")==0 )
    {
      scene2();
    }
    else if (m_params.scene().compare("3")==0 )
    {
      scene3();
    }
    else if (m_params.scene().compare("4")==0 )
    {
      scene4();
    }
    else if (m_params.scene().compare("5")==0 )
    {
      scene5();
    }
    else
    {
      std::cerr << "Application::run(): Could not initialize scene " << m_params.scene() << std::endl;
    }

    // Dump meta data so we know about setup and hardware etc.
    {
      std::string filename = m_params.output_path() + "/" + m_params.prefix() + "meta" + m_params.postfix() + ".cfg" ;

      std::ofstream file(filename.c_str());
      if (!file)
      {
        std::cerr << "Application::run(): Unable to open file " << filename.c_str() << std::endl;
        return;
      }

      file << "Experiment Parameters:" << std::endl;
      file << print_params(m_params);

      file << "Computational software:" << std::endl;
      file << print_versions() << std::endl;

      file << "Computational device:" << std::endl;
      file << print_device(m_params.device_number()) << std::endl;

      file.flush();
      file.close();
    }

    m_cur_time = 0.0;

    while( m_cur_time < m_params.total_time() )
    {
      ::run(m_simulator, m_params);

      m_cur_time = m_cur_time + m_params.time_step();

      SVGFile file_svg;
      SVGFile * svg = 0;

      std::stringstream filename;
      static int cnt = 0;
      filename << m_params.output_path() <<  "/" << m_params.prefix() << "image" << m_params.postfix() <<"_" << std::setw(4) << std::setfill('0') << ++cnt << ".svg";
      
      file_svg.open(filename.str(), 1.0, 1.0, 500);
      
      svg = &file_svg;

      draw(
           m_simulator
           , m_params.draw_grid()
           , m_params.draw_walls()
           , m_params.draw_liquid()
           , m_params.draw_particles()
           , m_params.draw_velocities()
           , m_params.draw_unilateral()
           , m_params.draw_solution()
           , m_params.draw_divergence()
           , svg);
    }

    std::cout << "Application::run(): Whole simulation is done: " << m_cur_time << " secs" << std::endl;

    if(m_params.profiling())
    {
      std::string filename = m_params.output_path() + "/" + m_params.prefix() + "profiling" + m_params.postfix() + ".m" ;

      std::ofstream file;

      file.open(filename.c_str());

      file << matlab_write_profiling() << std::endl;

      file.flush();
      file.close();

      std::cout << "Application::run() Done writing profile data..." << std::endl;
    }

  }

};


int main(int argc, char **argv)
{
  Application app;
  app.run(argc,argv);

	return 0;
}
