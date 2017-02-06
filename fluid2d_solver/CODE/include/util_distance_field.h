#ifndef UTIL_DISTANCE_FIELD_H
#define UTIL_DISTANCE_FIELD_H

#include <util_vec.h>

#include <cmath>

template<typename F>
class DistanceFieldExp
{
public:
  
  float operator()(Vec2f const & p) const 
  {
    F const & self = static_cast<F const &>( *this );
    return self(p);
  }

  operator F &()             { return static_cast<      F&>(*this); }
  operator F const &() const { return static_cast<const F&>(*this); }
  
};

class Plane
: public DistanceFieldExp<Plane>
{
public:
  
  float m_nx;
  float m_ny;
  float m_w;
  
public:
  
  Plane()
  : m_nx(1.0f)
  , m_ny(0.0f)
  , m_w(0.5)
  {
  }
  
  Plane(float const nx, float const ny,float const w)
  {
    using std::sqrt;
    
    m_nx = nx / sqrt(nx*nx + ny*ny);
    m_ny = ny / sqrt(nx*nx + ny*ny);
    m_w  = w;
  }
    
  float operator()(Vec2f const & p) const
  {
    float const & x = p[0];
    float const & y = p[1];
    
    return (x*m_nx + y*m_ny - m_w);
  }
  
};

class Box
	: public DistanceFieldExp<Box>
{
public:
  
  float m_min_x;
  float m_max_x;
  float m_min_y;
  float m_max_y;
  
public:
  
  Box()
  : m_min_x(0.1)
  , m_max_x(0.9)
  , m_min_y(0.1)
  , m_max_y(0.9)
  {
  }
  
  Box(float const min_x,float const max_x,float const min_y,float const max_y)
  : m_min_x(min_x)
  , m_max_x(max_x)
  , m_min_y(min_y)
  , m_max_y(max_y)
  {
  }
  
  virtual ~Box()
  {
  }
  
  float operator()(Vec2f const & p) const
  {
    using std::fabs;
    using std::sqrt;
    
    float const & x = p[0];
    float const & y = p[1];
    
    float const cx = (m_min_x+m_max_x)/2.0f;  // box center
    float const cy = (m_min_y+m_max_y)/2.0f;
    float const w  = (m_max_x-m_min_x)/2.0f;  // half widht and height
    float const h  = (m_max_y-m_min_y)/2.0f;
    
    float const dx = fabs(x-cx);       // point in positive ortant
    float const dy = fabs(y-cy);
    float const qx = (dx>w) ? w : dx;        // project onto box if point is outside
    float const qy = (dy>h) ? h : dy;
    
    float const dist2 = (qx-dx)*(qx-dx) + (qy-dy)*(qy-dy) ;  // squared distance to closest point on surface
    
    float const tx =  w-qx;  // if inside then distance to box wall
    float const ty =  h-qy;
    
    return dist2 > 0 ? - sqrt(dist2) : (tx<ty ? tx : ty);  
  }
  
};

class Circle
: public DistanceFieldExp<Circle>
{
public:
  
  float m_x;
  float m_y;
  float m_r;
  
public:
  
  Circle()
  : m_x(0.5)
  , m_y(0.5)
  , m_r(0.4)
  {
  }
  
  Circle(float const x,float const y,float const r)
  : m_x(x)
  , m_y(y)
  , m_r(r)
  {
  }
  
  virtual ~Circle()
  {
  }
  
  float operator()(Vec2f const & p) const
  {    
    return - (dist(p,Vec2f(m_x,m_y)) - m_r);
  }  
};

template<typename F1, typename F2>
class Intersection
: public DistanceFieldExp< Intersection<F1,F2> >
{
public:
  
  F1 const & m_A;
  F2 const & m_B;

public:
  
  Intersection( DistanceFieldExp<F1> const & A, DistanceFieldExp<F2> const & B)
  : m_A(A)
  , m_B(B)
  {
  }
  
  virtual ~Intersection()
  {
  }

public:
  
  float operator()(Vec2f const & p) const
  {    
    using std::min;
    
    return min( m_A(p), m_B(p) );
  }  
  
};

template<typename F>
class Invert
: public DistanceFieldExp< Invert<F> >
{
public:
  
  F const & m_A;
  
public:
  
  Invert(DistanceFieldExp<F> const & A)
  : m_A(A)
  {
  }
  
  virtual ~Invert()
  {
  }
  
public:
  
  float operator()(Vec2f const & p) const
  {    
    return - m_A(p);
  }  
  
};

template<typename F1, typename F2>
class Union
: public DistanceFieldExp< Union<F1,F2> >
{
public:
  
  F1 const & m_A;
  F2 const & m_B;
  
public:
  
  Union(DistanceFieldExp<F1> const & A, DistanceFieldExp<F2> const & B)
  : m_A(A)
  , m_B(B)
  {
  }
    
public:
  
  float operator()(Vec2f const & p) const
  {    
    using std::max;
    
    return max( m_A(p), m_B(p) );
  }  
  
};

class Const
: public DistanceFieldExp<Const>
{
public:
  
  float m_val;
  
public:
  
  Const()
  : m_val(0.0f)
  {
  }
  
  Const(float const val)
  : m_val(val)
  {
  }
  
  float operator()(Vec2f const & p) const
  {
    return m_val;
  }
  
};

template<typename F>
class Dilation
: public DistanceFieldExp< Dilation<F> >
{
public:
  
  F     const & m_A;
  Const         m_B;
  
public:
  
  Dilation(DistanceFieldExp<F> const & A, Const const & B)
  : m_A(A)
  , m_B(B)
  {
  }
  
  float operator()(Vec2f const & p) const
  {
    return m_A(p) + m_B(p);
  }
  
};

template <typename F1, typename F2>
inline Union<F1,F2> const operator+(DistanceFieldExp<F1> const& A, DistanceFieldExp<F2> const& B)
{
  return Union<F1,F2>(A,B);
}

template <typename F1, typename F2>
inline Intersection<F1,F2> const operator%(DistanceFieldExp<F1> const& A, DistanceFieldExp<F2> const& B)
{
  return Intersection<F1,F2>(A,B);
}

template <typename F>
inline Invert<F> const operator-(DistanceFieldExp<F> const& A)
{
  return Invert<F>(A);
}

template <typename F>
inline Dilation<F> const dilation(DistanceFieldExp<F> const& A, float const & val)
{
  return Dilation<F>(A,Const(val));
}

template <typename F>
inline Dilation<F> const dilation(float const & val,DistanceFieldExp<F> const& A)
{
  return A+val;
}


// UTIL_DISTANCE_FIELD_H
#endif
