#ifndef PROFILING_H
#define PROFILING_H

#include <util_timer.h>

#include <map>
#include <vector>
#include <cmath>
#include <sstream>
#include <cassert>

#define START_TIMER(NAME) Profiling::get_timer_monitor(NAME)->start()
#define STOP_TIMER(NAME) Profiling::get_timer_monitor(NAME)->stop()
#define RECORD(NAME,VALUE) Profiling::get_monitor(NAME)->record(VALUE)
#define NEW_UNIQUE_NAME(NAME) Profiling::generate_unique_name(NAME)
#define RECORD_VECTOR(NAME,VALUE) Profiling::get_vector_monitor(NAME)->record(VALUE)

/**
 * Profiling Class.
 * This class provides functionality for monitoring/handling multiple
 * timers and recording profiling data.
 */
class Profiling
{
public:

  class Monitor
  {
  protected:

    float       m_total;
    float       m_min;
    float       m_max;
    float       m_avg;
    size_t      m_N;
    std::vector<float> m_values;

  public:

    float get_avg() const    { return m_avg;    }
    float get_min() const    { return m_min;    }
    float get_max() const    { return m_max;    }
    float get_total() const  { return m_total;  }
    std::vector<float> const & get_values() const { return m_values; }
    size_t get_number_of_entries() const  { return m_N;  }

  public:

    Monitor()
    : m_total( 0.0f )
    , m_min( 0.0f )
    , m_max( 0.0f )
    , m_N( 0u )
    , m_values()
    {}

  public:

    void record( float const & value )
    {
      using std::min;
      using std::max;

      m_values.push_back( value );
      m_min   = min( m_min, value );
      m_max   = max( m_max, value );
      m_avg   =  (m_N*m_avg + value)/ (m_N+1);
      m_N     = m_N + 1;
      m_total = m_total + value;
    }
    void clear()
    {
      using std::min;
      using std::max;

      m_values.clear();
      m_min   = 0.0f;
      m_max   = 0.0f;
      m_avg   = 0.0f;
      m_N     = 0u;
      m_total = 0.0f;
    }
  };

  class TimerMonitor
  : public Monitor
  {
  protected:

    Timer m_timer;

  public:

    TimerMonitor()
    : Monitor()
    , m_timer()
    {}

  public:

    void start()
    {
      m_timer.start();
    }

    void stop()
    {
      m_timer.stop();
      float msecs = m_timer();
      Monitor::record( msecs );
    }

  };

  class VectorMonitor
  {
  protected:

    std::vector< std::vector<float>  > m_values;

  public:

    std::vector<float> const & get_values(size_t const & idx) const { return m_values[idx]; }
    size_t get_number_of_entries() const  { return m_values.size();  }

  public:

    VectorMonitor()
    : m_values()
    {}

  public:

    void record( std::vector<float> & value )
    {
      m_values.push_back( value );
    }

    void clear()
    {
      using std::min;
      using std::max;

      m_values.clear();
    }
  };

  typedef std::map<std::string, Monitor >        monitors_container;
  typedef std::map<std::string, TimerMonitor>    timer_monitors_container;
  typedef std::map<std::string, VectorMonitor >  vector_monitors_container;
  typedef std::map<std::string, size_t >         counter_container;

private:



  static counter_container & counters()
  {
    static counter_container data;

    return data;
  }

public:

  static std::string generate_unique_name(std::string const & name)
  {
    size_t const value = counters()[name]++;

    std::ostringstream s;

    s << name << "_" << value;

    return s.str();
  }

  static monitors_container * get_monitors_instance()
  {
    static monitors_container  monitors;
    return &monitors;
  }

  static timer_monitors_container * get_timer_monitors_instance()
  {
    static timer_monitors_container  monitors;
    return &monitors;
  }

  static vector_monitors_container * get_vector_monitors_instance()
  {
    static vector_monitors_container  monitors;
    return &monitors;
  }

  static Monitor * get_monitor(std::string const & name)
  {
    return &(*Profiling::get_monitors_instance())[name];
  }

  static TimerMonitor * get_timer_monitor(std::string const & name)
  {
    return &(*Profiling::get_timer_monitors_instance())[name];
  }

  static VectorMonitor * get_vector_monitor(std::string const & name)
  {
    return &(*Profiling::get_vector_monitors_instance())[name];
  }

};


// PROFILING_H
#endif
