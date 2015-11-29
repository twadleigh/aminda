#include "vec.hpp"

extern "C"
void map_to_plane(
  float const *homog, int nc,
  int snx, int sny, float const *src,
  int tnx, int tny, float *tgt)
{
  mat3f const &hom = *reinterpret_cast<mat3f const *>(homog);

  int const sdc = 1;
  int const sdx = sdc*nc;
  int const sdy = sdx*snx;
  int const tdc = 1;
  int const tdx = tdc*nc;
  int const tdy = tdx*tnx;

#pragma omp parallel for
  for (int j = 0; j < tny; ++j) {
    for (int i = 0; i < tnx; ++i) {
      float *_tgt = tgt + i*tdx + j*tdy;

      vec3f p = mvmult(hom, vec3f(float(i)+0.5,float(j)+0.5,1.0));
      float x = p[0]/p[2];
      float y = p[1]/p[2];
      if (x > 0.0f && x < snx && y > 0.0f && y < sny) {
        int xlo(x), ylo(y);
        float bx = x-float(xlo); float ax = 1.0f-bx;
        float by = y-float(ylo); float ay = 1.0f-by;
        float aa = ax*ay; float ab = ax*by;
        float ba = bx*ay; float bb = bx*by;
        float const *saa = src + xlo*sdx + ylo*sdy;
        float const *sab = saa + sdy;
        float const *sba = saa + sdx;
        float const *sbb = sba + sdy;
        for(int c = 0; c < nc; ++c) {
          _tgt[c] = aa*saa[c]+ab*sab[c]+ba*sba[c]+bb*sbb[c];
        }
      } else {
        for (int c = 0; c < nc; ++c) _tgt[c] = 0.0f;
      }
    }
  }
}

extern "C"
void mean_and_deviation(
  int w, int nc, int nx, int ny,
  float const *image, float *mean)
{
  int const dc = 1;
  int const dx = dc*nc;
  int const dy = dx*nx;

  float const invScale = 1.0f/float((2*w+1)*(2*w+1));

#pragma omp parallel for
  for (int j = w; j < ny-w; ++j) {
    for (int i = w; i < nx-w; ++i) {
      float *_mean = mean + i*dx + j*dy;

      // accumulate mean
      for (int jj = j-w; jj <= j+w; ++jj) {
        float const *_image_y = image + jj*dy;
        for (int ii = i-w; ii <= i+w; ++ii) {
          float const *_image = _image_y + ii*dx;
          for (int c = 1; c < nc; ++c) {
            _mean[c] += _image[c];
          }
        }
      }
      for (int c = 1; c < nc; ++c) _mean[c] *= invScale;

      // accumulate deviation
      for (int jj = j-w; jj <= j+w; ++jj) {
        float const *_image_y = image + jj*dy;
        for (int ii = i-w; ii <= i+w; ++ii) {
          float const *_image = _image_y + ii*dx;
          for (int c = 1; c < nc; ++c) {
            float diff = _image[c] - _mean[c];
            _mean[0] += diff*diff;
          }
        }
      }
      _mean[0] = sqrt(_mean[0]);
    }
  }
}

extern "C"
void normalized_cross_correlation(
  int w, int nc, int nx, int ny,
  float const *im1, float const *mn1, float const *im2, float const *mn2,
  float *nxcorr)
{
  int const dx = nc;
  int const dy = dx*nx;

  for (int j = w; j < ny-w; ++j) {
    for (int i = w; i < nx-w; ++i) {
      int ix = i+j*nx;
      int _ix = ix*nc;

      float const *_mn1 = mn1+_ix;
      float const *_mn2 = mn2+_ix;
      float const *_im1 = im1+_ix;
      float const *_im2 = im2+_ix;

      if(_im1[0] == 0.0f || _im2[0] == 0.0f) {
        nxcorr[ix] = -1.0f;
      } else {
        nxcorr[ix] = 0.0f;
        for (int jj = j-w; jj <= j+w; ++jj) {
          int ojj = jj*dy;
          float const *_im1jj = im1 + ojj;
          float const *_im2jj = im2 + ojj;
          for (int ii = i-w; ii <= i+w; ++ii) {
            int oii = ii*dx;
            float const *_i1 = _im1jj + oii;
            float const *_i2 = _im2jj + oii;
            for (int c = 1; c < nc; ++c) {
              nxcorr[ix] += (_i1[c] - _mn1[c]) * (_i2[c] - _mn2[c]);
            }
          }
        }
        nxcorr[ix] /= mn1[0] * mn2[0];
      }
    }
  }
}

extern "C"
void fill(int nx, int ny, float *a, float x) {
  int const dx = 1;
  int const dy = nx;

#pragma omp parallel for
  for (int j = 0; j < ny; ++j) {
    float *_a = a + j*dy;
    for (int i = 0; i < nx; ++i) {
      _a[i] = x;
    }
  }
}

extern "C"
void maximum(int nx, int ny, float const *a, float const *b, float *mx) {
  int const dx = 1;
  int const dy = nx;

#pragma omp parallel for
  for (int j = 0; j < ny; ++j) {
    int offset = j*dy;
    float const *_a = a + offset;
    float const *_b = b + offset;
    float *_mx = mx + offset;
    for (int i = 0; i < nx; ++i) {
      float aa = _a[i];
      float bb = _b[i];
      _mx[i] = aa > bb ? aa : bb;
    }
  }
}

extern "C"
void update_depth(int nx, int ny, float const *nxc, float *mxnxc, float *zs, float z) {
  int const dx = 1;
  int const dy = nx;

#pragma omp parallel for
  for (int j = 0; j < ny; ++j) {
    int offset = j*dy;
    float const *_nxc = nxc + offset;
    float *_mxnxc = mxnxc + offset;
    float *_zs = zs + offset;
    for (int i = 0; i < nx; ++i) {
      float n = _nxc[i];
      if(n > _mxnxc[i]) {
        _mxnxc[i] = n;
        _zs[i] = z;
      }
    }
  }
}
