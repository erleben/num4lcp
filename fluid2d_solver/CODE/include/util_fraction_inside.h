#ifndef UTIL_FRACTION_INSIDE_H
#define UTIL_FRACTION_INSIDE_H

/**
 * Given two signed distance values, determine what fraction of a connecting
 * segment is "inside"
 */
template<typename T>
inline T fraction_inside(T const & phi_left, T const & phi_right)
{
  if(phi_left < 0 && phi_right < 0)
    return 1;

  if (phi_left < 0 && phi_right >= 0)
    return phi_left / (phi_left - phi_right);

  if(phi_left >= 0 && phi_right < 0)
    return phi_right / (phi_right - phi_left);

  return 0;
}

// UTIL_FRACTION_INSIDE_H
#endif
