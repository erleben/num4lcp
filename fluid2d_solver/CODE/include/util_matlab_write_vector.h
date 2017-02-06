#ifndef UTIL_MATLAB_WRITE_VECTOR_H
#define UTIL_MATLAB_WRITE_VECTOR_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

template<typename T>
inline std::string matlab_write_vector(std::vector<T> const & values, size_t const & size, bool const & add_one = false)
{
  std::stringstream output;

  output << "[";

  if(add_one)
  {
    for(size_t i = 0u; i< size; ++i)
    {
      output <<  values[i] + T(1) << " ";
    }
  }
  else
  {
    for(size_t i = 0u; i< size; ++i)
    {
      output <<  values[i] << " ";
    }
  }

  output << "]";

  output.flush();

  return output.str();
}

template<typename T>
inline std::string matlab_write_vector(std::vector<T> const & values )
{
  return matlab_write_vector(values, values.size(), false);
}


template<typename T>
inline std::string matlab_write_vector(std::string const & name
                                       , std::vector<T> const & values)
{
  std::stringstream output;

  output << name << " = [";

  for(size_t i = 0u; i< values.size(); ++i)
  {
    output << values[i] << " ";
  }
  output << "];";

  output.flush();

  return output.str();
}

// UTIL_MATLAB_WRITE_VECTOR_H
#endif
