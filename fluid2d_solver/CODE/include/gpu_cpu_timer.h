#ifndef GPU_CPU_TIMER_H
#define GPU_CPU_TIMER_H

#include <time.h>

/**
* A CPU Timer.
* This class is intended to make it more easy to use CPU timers.
* It measures the time duration between a start and stop call in milliseconds. Observe that everything is measured inbetween the start and stop calls, like the wall clock.
*/
class CPUTimer
{
protected:

	clock_t         m_start;    ///<  A data structure to identify the start event.
	clock_t         m_end;      ///<  A data structure to identify the start event.

public:

	void start() { m_start = clock(); }

	void stop()  { m_end = clock(); }

	/**
	* Get Time.
	*
	* @return     The time between start and stop in milliseconds.
	*/
	float operator()()const
	{
	  return ( m_end - m_start ) / ((float) CLOCKS_PER_SEC) * 1000.0f;
	}
};

// GPU_CPU_TIMER_H
#endif 
