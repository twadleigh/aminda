#pragma once

namespace aminda {

  template<size_t N, typename F>
  class vec {
  public:
    const F& operator[](int ix) const { return _data[ix]; }
    F& operator[](int ix) { return _data[ix]; }
  private:
    F _data[N];
  };

}
