#ifndef UTIL_SVG_FILE_H
#define UTIL_SVG_FILE_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>

#include <cassert>

struct SVGStyle
{
public:

  std::string m_stroke;
  std::string m_fill;
  double      m_stroke_width;

public:

  SVGStyle()
  : m_stroke("black")
  , m_fill("black")
  , m_stroke_width(0.5)
  {
  }

  SVGStyle(std::string const & stroke, std::string const & fill)
  : m_stroke(stroke)
  , m_fill(fill)
  , m_stroke_width(0.5)
  {
  }


  static std::string make_hex_color(double const & red, double const & green, double const & blue)
  {
    assert( red   <= 1.0 || !"illegal red color value"  );
    assert( red   >= 0.0 || !"illegal red color value"  );
    assert( green <= 1.0 || !"illegal green color value");
    assert( green >= 0.0 || !"illegal green color value");
    assert( blue  <= 1.0 || !"illegal blue color value" );
    assert( blue  >= 0.0 || !"illegal blue color value" );

    unsigned int const int_red   =   0xFF*red;
    unsigned int const int_green =   0xFF*green;
    unsigned int const int_blue  =   0xFF*blue;

    std::ostringstream ss;
    ss << "#";
    ss.width(2);
    ss.fill('0');
    ss << std::hex << int_red;
    ss.width(2);
    ss.fill('0');
    ss << std::hex << int_green;
    ss.width(2);
    ss.fill('0');
    ss << std::hex << int_blue;
    return ss.str();
  }

};

/**
 * Have a look at http://www.svgbasics.com
 */
class SVGFile
{
protected:

  std::ofstream m_stream;

  double m_dx;
  double m_dy;
  double m_scale;
  unsigned int m_arrow_count;

  SVGStyle m_style;

public:

  SVGFile()
  : m_stream()
  , m_dx(300.0)
  , m_dy(300.0)
  , m_scale(2.0)
  , m_style()
  {
  }

  ~SVGFile()
  {
    if(m_stream.is_open() )
      this->close();
  }

protected:

  double dx(double x) { return x * m_scale; }
  double dy(double y) { return (m_dy - y) * m_scale; } // svg has inverted y axis

  std::string make_style_attr(bool stroke = true)
  {
    std::ostringstream ss;
    ss << "style=\"";
    if (stroke)
    {
      ss
      //<< "stroke-linecap:round;" // Does not work for lines
      << "stroke-width:" << m_style.m_stroke_width << ";"
      << "stroke:" << m_style.m_stroke << ";";
    }
    ss << "fill:" << m_style.m_fill << ";\"";
    return ss.str();
  }

  void close()
  {
    assert(m_stream.is_open() || !"internal error");

    m_stream << "</svg>" << std::flush;
    m_stream.flush();
    m_stream.close();
  }

public:

  void open(std::string const & filename, double const & width, double const & height, double const & scale = 1.0)
  {
    assert( width > 0.0 || !"width must be positive");
    assert( height > 0.0 || !"height must be positive");
    assert( height > 0.0 || !"scale must be positive");

    m_stream.open( filename.c_str() );

    assert(m_stream.is_open() || !"internal error");

    //--- Set bounding box information
    // y axis is inverted in svg, so translate to negative y (and then it will be reversed when drawing)
    this->m_dx    = width;
    this->m_dy    = height;
    this->m_scale = scale;
    this->m_arrow_count = 0;

    m_stream << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
    m_stream << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
    << "width=\"" << m_scale*width << "\" "
    << "height=\"" << m_scale*height << "\" "
    << "version=\"1.1\"> "
    << std::endl;
  }


  void draw_arrow2d(double const & x1, double const & y1,
                    double const & x2, double const & y2,
                    double const xscale = 1.5, double const yscale = 1.0)
  {
    assert(m_stream.is_open() || !"internal error");

    double xmax = xscale*10.0;
    double ymax = yscale*10.0;
    double refX = xmax/10;
    double refY = ymax/2;

    m_stream << "<defs>"
        << "<marker id=\"my_arrow_" << this->m_arrow_count << "\" "
        << "viewBox=\"0 0 " << xmax << " " << ymax << "\" "
        << "refX=\"" << refX << "\" refY=\"" << refY << "\" "
        << "markerUnits=\"strokeWidth\" orient=\"auto\" "
        << "markerWidth=\"" << refX*3 << "\" "                  
        << "markerHeight=\"" << refX*3 << "\">"
        << "<polyline points=\"0,0 " << xmax << "," << refY << " 0," << ymax << " " << refX << ","<< refY << "\" "
        << "fill=\"" << this->m_style.m_fill << "\" />"
        << "</marker>"
        << "</defs>"
        << "<line "
        << "x1=\"" << dx(x1) << "\" "
        << "y1=\"" << dy(y1) << "\" "
        << "x2=\"" << dx(x2) << "\" "
        << "y2=\"" << dy(y2) << "\" "
        << "stroke=\"" << this->m_style.m_stroke << "\" "
        << "stroke-width=\"" << this->m_style.m_stroke_width << "\" "
        << "marker-end=\"url(#my_arrow_" << this->m_arrow_count << ")\" "
        << " />" << std::endl;

    this->m_arrow_count++;

  }

  void draw_line(double const & x1, double const & y1, double const & x2, double const & y2)
  {
    assert(m_stream.is_open() || !"internal error");

    m_stream << "<line "
        << "x1=\"" << dx(x1) << "\" "
        << "y1=\"" << dy(y1) << "\" "
        << "x2=\"" << dx(x2) << "\" "
        << "y2=\"" << dy(y2) << "\" "
        << make_style_attr() << " />" << std::endl;
  }

  void draw_circle(double const & x, double const & y, double const & r)
  {
    assert(m_stream.is_open() || !"internal error");

    m_stream << "<circle "
    << "cx=\"" << dx(x) << "\" "
    << "cy=\"" << dy(y) << "\" "
    << "r=\"" << this->m_scale*r << "\" "
    << make_style_attr() << " />" << std::endl;
  }

  void draw_rect(double const & x, double const & y, double const & w,double const & h)
  {
    assert(m_stream.is_open() || !"internal error");

    m_stream << "<rect "
    << "x=\"" << dx(x - (w/2.0)) << "\" "
    << "y=\"" << dy(y - (h/2.0)) << "\" "
    << "width=\"" << this->m_scale*w << "\" "
    << "height=\"" << this->m_scale*h << "\" "
    << make_style_attr() << " />" << std::endl;
  }

  void set_color(double const & red, double const & green, double const & blue)
  {
    set_color(   SVGStyle::make_hex_color(red,green,blue) );
  }

  void set_color(std::string color)
  {
      m_style.m_fill = color;
      m_style.m_stroke = color;
  }

  void set_fill(double const & red, double const & green, double const & blue)
  {
    std::ostringstream color;
    color << "rgb(" << red << "," << green << "," << blue << ")";

    set_fill(color.str());
  }

  void set_fill(std::string color)
  {
    m_style.m_fill = color;
  }

  void set_stroke(double const & red, double const & green, double const & blue)
  {
    std::ostringstream color;
    color << "rgb(" << red << "," << green << "," << blue << ")";

    set_stroke(color.str());
  }

  void set_stroke(std::string color)
  {
    m_style.m_stroke = color;
  }

  void set_stroke_width(double const & width)
  {
    m_style.m_stroke_width  = width;
  }
  
};


// UTIL_SVG_FILE_H
#endif
