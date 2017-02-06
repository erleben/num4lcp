#ifndef GPU_CUDA_TIMER_H
#define GPU_CUDA_TIMER_H

#include <cuda_runtime.h>    // needed for CUDA C++ runtime api

#include <gpu_safe_call.h>    // needed for CUDA C++ runtime api


/**
 * A CUDA Timer.
 * This class wraps cuda runtime functionality for measuring the time spend between a start and stop event.
 * The intention is that it makes it easier to write code for measuring time spendure in kernel wrappers.
 */
class CUDATimer
  {
  protected:
    
    cudaEvent_t m_start;   ///< Cuda runtime data structure used to identify the start event.
    cudaEvent_t m_stop;    ///< Cuda runtime data structure used to identify the stop event.
    float       m_time;    ///< Variable used to hold the time between the last start and stop event in units of milliseconds.

  public:
    
    CUDATimer()
    : m_start()
    , m_stop()
    , m_time(0.0f)
    {
    }
    
    ~CUDATimer()
    {
    }
    
  public:
    
    void start()
    {
      SAFE_CALL( cudaEventCreate( &m_start ) ); 
      SAFE_CALL( cudaEventCreate( &m_stop  ) ); 
      SAFE_CALL( cudaEventRecord( m_start, 0 ) );
    }
    
    void stop()
    {
      SAFE_CALL( cudaEventRecord( m_stop, 0 ) );
      SAFE_CALL( cudaEventSynchronize( m_stop ) ); 
      SAFE_CALL( cudaEventElapsedTime( &m_time, m_start, m_stop ) ); 
      SAFE_CALL( cudaEventDestroy( m_start ) ); 
      SAFE_CALL( cudaEventDestroy( m_stop )  ); 
    }
    
    float milliseconds() const { return m_time; }
    
    /**
     * Get Time.
     *
     * @return       The time between last start and stop events in milliseconds.
     */
    float operator()() const {  return milliseconds();  }
      
  };

// GPU_CUDA_TIMER_H
#endif 
