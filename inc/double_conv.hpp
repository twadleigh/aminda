#pragma once

#include <cstdint>

namespace aminda {

  template<typename T>
  inline T fromDouble(double d) { return static_cast<T>(d); }

  template<typename T>
  inline double toDouble(T t) { return static_cast<double>(t); }

  template<>
  inline double fromDouble<double>(double d) { return d; }

  template<>
  inline double toDouble<double>(double d) { return d; }

  template<>
  inline uint8_t fromDouble<uint8_t>(double d) { 
    const double scale = static_cast<double>(UINT8_MAX);
    return static_cast<uint8_t>(scale*d); 
  }

  template<>
  inline double toDouble<uint8_t>(uint8_t t) { 
    const double scale = 1.0/static_cast<double>(UINT8_MAX);
    return static_cast<double>(t)*scale;
  }

  template<>
  inline uint16_t fromDouble<uint16_t>(double d) { 
    const double scale = static_cast<double>(UINT16_MAX);
    return static_cast<uint16_t>(scale*d); 
  }

  template<>
  inline double toDouble<uint16_t>(uint16_t t) { 
    const double scale = 1.0/static_cast<uint16_t>(UINT16_MAX);
    return static_cast<double>(t)*scale;
  }

}
