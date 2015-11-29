#pragma once

#include <cstring>
#include <cmath>

template<size_t N, typename F>
class vec {
public:
  vec(const F& f = F()) {
    for (size_t i = 0; i < N; ++i)
      _d[i] = f;
  }

  template<size_t M, typename U>
  vec(const vec<M, U>& v, const F& f = F()) {
    for (size_t i = 0; i < M && i < N; ++i)
      _d[i] = static_cast<F>(v[i]);
    for (size_t i = M; i < N; ++i)
      _d[i] = f;
  }

  template<size_t M, typename U>
  vec(const F& f, const vec<M, U>& v) {
    static const size_t DIFF = N > M ? N - M : 0;
    for (size_t i = 0; i < DIFF; ++i)
      _d[i] = f;
    for (size_t i = DIFF; i < N; ++i)
      _d[i] = static_cast<F>(v[i - DIFF]);
  }

  vec(const F& v0, const F& v1) {
    _d[0] = v0;
    _d[1] = v1;
    for (size_t i = 2; i < N; ++i)
      _d[i] = F();
  }

  vec(const F& v0, const F& v1, const F& v2) {
    _d[0] = v0;
    _d[1] = v1;
    _d[2] = v2;
    for (size_t i = 3; i < N; ++i)
      _d[i] = F();
  }

  vec(const F& v0, const F& v1, const F& v2, const F& v3) {
    _d[0] = v0;
    _d[1] = v1;
    _d[2] = v2;
    _d[3] = v3;
    for (size_t i = 4; i < N; ++i)
      _d[i] = F();
  }

  vec(const F& v0, const F& v1, const F& v2, const F& v3, const F& v4) {
    _d[0] = v0;
    _d[1] = v1;
    _d[2] = v2;
    _d[3] = v3;
    _d[4] = v4;
    for (size_t i = 5; i < N; ++i)
      _d[i] = F();
  }

  vec(const F& v0, const F& v1, const F& v2, const F& v3, const F& v4,
      const F& v5) {
    _d[0] = v0;
    _d[1] = v1;
    _d[2] = v2;
    _d[3] = v3;
    _d[4] = v4;
    _d[5] = v5;
    for (size_t i = 6; i < N; ++i)
      _d[i] = F();
  }

  vec(const F& v0, const F& v1, const F& v2, const F& v3, const F& v4,
      const F& v5, const F& v6) {
    _d[0] = v0;
    _d[1] = v1;
    _d[2] = v2;
    _d[3] = v3;
    _d[4] = v4;
    _d[5] = v5;
    _d[6] = v6;
    for (size_t i = 7; i < N; ++i)
      _d[i] = F();
  }

  vec(const F& v0, const F& v1, const F& v2, const F& v3, const F& v4,
      const F& v5, const F& v6, const F& v7) {
    _d[0] = v0;
    _d[1] = v1;
    _d[2] = v2;
    _d[3] = v3;
    _d[4] = v4;
    _d[5] = v5;
    _d[6] = v6;
    _d[7] = v7;
    for (size_t i = 8; i < N; ++i)
      _d[i] = F();
  }

  template<typename Ix>
  const F& operator[](Ix ix) const {
    return _d[ix];
  }

  template<typename Ix>
  F& operator[](Ix ix) {
    return _d[ix];
  }

  vec<N, F>& operator+=(const vec<N, F>& v) {
    for (size_t i = 0; i < N; ++i)
      _d[i] += v._d[i];
    return *this;
  }

  vec<N, F>& operator-=(const vec<N, F>& v) {
    for (size_t i = 0; i < N; ++i)
      _d[i] -= v._d[i];
    return *this;
  }

  vec<N, F>& operator*=(const vec<N, F>& v) {
    for (size_t i = 0; i < N; ++i)
      _d[i] *= v._d[i];
    return *this;
  }

  vec<N, F>& operator*=(const F& f) {
    for (size_t i = 0; i < N; ++i)
      _d[i] *= f;
    return *this;
  }

  vec<N, F>& operator/=(const vec<N, F>& v) {
    for (size_t i = 0; i < N; ++i)
      _d[i] /= v._d[i];
    return *this;
  }

  vec<N, F>& operator/=(const F& f) {
    for (size_t i = 0; i < N; ++i)
      _d[i] /= f;
    return *this;
  }

	vec<N, F> operator-() const {
		vec<N, F> v;
		for (size_t i = 0; i < N; ++i)
			v._d[i] = -_d[i];
		return v;
	}

  vec<N, F> operator+(const vec<N, F>& v) const {
    return vec<N, F>(*this) += v;
  }

  vec<N, F> operator-(const vec<N, F>& v) const {
    return vec<N, F>(*this) -= v;
  }

  vec<N, F> operator*(const F& f) const {
    return vec<N, F>(*this) *= f;
  }

  vec<N, F> operator*(const vec<N, F>& v) const {
    return vec<N, F>(*this) *= v;
  }

  vec<N, F> operator/(const F& f) const {
    return vec<N, F>(*this) /= f;
  }

  vec<N, F> operator/(const vec<N, F>& v) const {
    return vec<N, F>(*this) /= v;
  }

private:
  F _d[N];
};

template<size_t N, typename F>
inline vec<N, F> operator+(const F& f, const vec<N, F>& v) {
  return v + f;
}

template<size_t N, typename F>
inline vec<N, F> operator-(const F& f, const vec<N, F>& v) {
  return vec<N, F>(f) - v;
}

template<size_t N, typename F>
inline vec<N, F> operator*(const F& f, const vec<N, F>& v) {
  return v * f;
}

template<size_t N, typename F>
inline vec<N, F> operator/(const F& f, const vec<N, F>& v) {
  return vec<N, F>(f) / v;
}

template<typename F>
inline F cross(const vec<2, F>& v1, const vec<2, F>& v2) {
  return v1[0] * v2[1] - v1[1] * v2[0];
}

template<typename F>
inline vec<3, F> cross(const vec<3, F>& v1, const vec<3, F>& v2) {
  return vec<3, F>(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2],
      v1[0] * v2[1] - v1[1] * v2[0]);
}

template<size_t N, typename F>
inline F inner(const vec<N, F>& v1, const vec<N, F>& v2) {
  F f = F();
  for (size_t i = 0; i < N; ++i)
    f += v1[i] * v2[i];
  return f;
}

template<size_t N, typename F>
inline F length2(const vec<N, F>& v) {
  return inner(v, v);
}

template<size_t N, typename F>
inline F length(const vec<N, F>& v) {
  return sqrt(length2(v));
}

template<size_t N, typename F>
inline F cosine(const vec<N, F>& v1, const vec<N, F>& v2) {
  return inner(v1, v2) / (length(v1) * length(v2));
}

template<size_t N, typename F>
inline F angle(const vec<N, F>& v1, const vec<N, F>& v2) {
  return acos(cosine(v1, v2));
}

template<size_t N, size_t M, typename F>
vec<N, vec<M, F> > inline outer(const vec<N, F>& v1, const vec<M, F>& v2) {
  vec<N, vec<M, F> > out;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      out[i][j] = v1[i] * v2[j];
    }
  }
  return out;
}

template<size_t N, typename F>
vec<N, vec<N, F> > inline diag(const vec<N, F>& diag) {
  vec<N, vec<N, F> > d;
  for (size_t i = 0; i < N; ++i)
    d[i][i] = diag[i];
  return d;
}

template<size_t N, typename F>
vec<N, F> inline diag(const vec<N, vec<N, F> >& mat) {
  vec<N, F> d;
  for (size_t i = 0; i < N; ++i)
    d[i] = mat[i][i];
  return d;
}

template<size_t N, typename F>
vec<N + 1, vec<N + 1, F> > inline homog(
    const vec<N, vec<N + 1, F> >& mat) {
  vec<N + 1, vec<N + 1, F> > m(mat, vec<N + 1, F>());
  m[N][N] = F(1.0);
  return m;
}

template<size_t N, typename F>
vec<N, vec<N, F> > inline id() {
  return diag(vec<N, F>(static_cast<F>(1.0)));
}

template<size_t NR, size_t NC, typename F>
inline vec<NC, vec<NR, F> > tr(const vec<NR, vec<NC, F> >& mat) {
  vec<NC, vec<NR, F> > t;
  for (size_t i = 0; i < NR; ++i) {
    for (size_t j = 0; j < NC; ++j) {
      t[j][i] = mat[i][j];
    }
  }
  return t;
}

template<size_t NR, size_t NC, typename F>
inline vec<NR, F> mvmult(const vec<NR, vec<NC, F> >& m,
    const vec<NC, F>& v) {
  vec<NR, F> p;
  for (size_t i = 0; i < NR; ++i)
    p[i] = inner(m[i], v);
  return p;
}

template<size_t NR, size_t N, size_t NC, typename F>
inline vec<NR, vec<NC, F> > mmmult(const vec<NR, vec<N, F> >& m1,
    const vec<N, vec<NC, F> >& m2) {
  vec<NR, vec<NC, F> > prod;
  for (size_t i = 0; i < NR; ++i) {
    for (size_t j = 0; j < NC; ++j) {
      prod[i][j] = F();
      for (size_t k = 0; k < N; ++k) {
        prod[i][j] += m1[i][k] * m2[k][j];
      }
    }
  }
  return prod;
}

template<size_t N, typename F>
inline F det(const vec<N, vec<N, F> >& m) {
  F d = F(0.0);
  size_t i;
  F sgn;
  for (i = 0, sgn = F(1.0); i < N; ++i, sgn *= F(-1.0)) {
    // select sub array
    vec<N, vec<N, F> > sub;
    for (size_t j = 1; j < N; ++j) {
      size_t k, l;
      for (k = 0, l = 0; k < N; ++k, ++l) {
        if (k == i)
          ++k;
        sub[j][l] = m[j][k];
      }
    }
    // accumulate its determinant
    d += sgn * m[0][i] * det(sub);
  }
  return d;
}

// base case of the recursion
template<typename F>
inline F det(const vec<1, vec<1, F> >& m) {
  return m[0][0];
}

template<size_t N, typename F>
inline vec<N, vec<N, F> > minv(const vec<N, vec<N, F> >& m) {
  F mult;
  vec<N, vec<N, F> > mat = m;
  vec<N, vec<N, F> > inv = id<N, F>();

  // make upper triangular with 1s on the diagonal
  for (size_t i = 0; i < N; ++i) {
    mult = mat[i][i];
    mat[i] /= mult;
    inv[i] /= mult;
    for (size_t j = i + 1; j < N; ++j) {
      mult = mat[j][i];
      mat[j] -= mat[i] * mult;
      inv[j] -= inv[i] * mult;
    }
  }

  // make the identity
  for (size_t i = N - 1; i > 0; --i) {
    for (size_t j = 0; j < i; ++j) {
      F mult = mat[j][i];
      mat[j] -= mat[i] * mult;
      inv[j] -= inv[i] * mult;
    }
  }

  return inv;
}

typedef vec<2,int> vec2i;
typedef vec<3,int> vec3i;
typedef vec<4,int> vec4i;
typedef vec<2,size_t> vec2sz;
typedef vec<3,size_t> vec3sz;
typedef vec<4,size_t> vec4sz;

typedef vec<2,float> vec2f;
typedef vec<3,float> vec3f;
typedef vec<4,float> vec4f;
typedef vec<2,double> vec2d;
typedef vec<3,double> vec3d;
typedef vec<4,double> vec4d;

typedef vec<2,vec2f> mat22f;
typedef vec<3,vec2f> mat32f;
typedef vec<4,vec2f> mat42f;
typedef vec<2,vec3f> mat23f;
typedef vec<3,vec3f> mat33f;
typedef vec<4,vec3f> mat43f;
typedef vec<2,vec4f> mat24f;
typedef vec<3,vec4f> mat34f;
typedef vec<4,vec4f> mat44f;
typedef vec<2,vec2d> mat22d;
typedef vec<3,vec2d> mat32d;
typedef vec<4,vec2d> mat42d;
typedef vec<2,vec3d> mat23d;
typedef vec<3,vec3d> mat33d;
typedef vec<4,vec3d> mat43d;
typedef vec<2,vec4d> mat24d;
typedef vec<3,vec4d> mat34d;
typedef vec<4,vec4d> mat44d;

typedef mat22f mat2f;
typedef mat33f mat3f;
typedef mat44f mat4f;
typedef mat22d mat2d;
typedef mat33d mat3d;
typedef mat44d mat4d;
