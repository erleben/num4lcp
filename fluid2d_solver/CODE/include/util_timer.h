#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

#include <time.h>   // needed for clock_t type

/**
 * A CPU Timer.
 * This class is intended to make it more easy to use CPU timers.
 * It measures the time duration between a start and stop call in
 * milliseconds. Observe that everything is measured inbetween the
 * start and stop calls, like the wall clock.
 */
class Timer
{
protected:

  clock_t         m_start;    ///<  A data structure to identify the start event.
  clock_t         m_end;      ///<  A data structure to identify the start event.

public:

  void start()
  {
    m_start = m_end = clock();
    /*gettimeofday(&m_start, &m_tz);*/
  }

  void stop()
  {
    m_end = clock();
    /*gettimeofday(&m_end, &m_tz);*/
  }

  /**
   * Get Time.
   *
   * @return     The time between start and stop in milliseconds.
   */
  float operator()()const
  {
    // The gettimeofday() function has a resolution of microseconds -- we convert into milliseconds
    //return (m_end.tv_sec - m_start.tv_sec)*1000.0f + (m_end.tv_usec -  m_start.tv_usec)/1000.0f;

    // 2010-06-09 Kenny code review: why do we need the float cast? I thought an int multiplied by a
    //                               float would be type cast into a float automatically by the compiler?
    return ( m_end - m_start ) / ((float) CLOCKS_PER_SEC) * 1000.0f;
  }
  
};

// UTIL_TIMER_H
#endif
