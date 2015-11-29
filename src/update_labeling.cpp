#include "half.hpp"

using namespace half_float;

extern "C"
void update_labeling(
  int nits,
  int nx, int ny, int nz,
  half const *rho, half const *sigma,
  half *phi, half *d, half *b)
{
  int const dx = 1;
  int const dy = nx;
  int const dz = nx*ny;
  int const _dx = 3*dx;
  int const _dy = 3*dy;
  int const _dz = 3*dz;

#pragma omp parallel
  {
    for(int it = 0; it < nits; ++it) {

      // update d
#pragma omp for
      for (int k = 0; k < nz-1; ++k) {
        for (int j = 0; j < ny-1; ++j) {
          for (int i = 0; i < nx-1; ++i) {
            int ix = i*dx + j*dy + k*dz;
            int _ix = 3*ix;
            half *_d = d + _ix;
            half *_b = b + _ix;
            half const *_rho = rho + ix;
            half *_phi = phi + ix;

            _d[0] = _phi[dx] - _phi[0] + _b[0];
            _d[1] = _phi[dy] - _phi[0] + _b[1];
            _d[2] = _phi[dz] - _phi[0] + _b[2];

            half len = sqrt(_d[0]*_d[0] + _d[1]*_d[1] + _d[2]*_d[2]);

            if(len <= _rho[0]) {
              _d[0] = _d[1] = _d[2] = half(0.0);
            } else {
              half s = half(1.0) - rho[0]/len;
              _d[0] *= s;
              _d[1] *= s;
              _d[2] *= s;
            }
          }
        }
      } // update d

      // update phi (even)
#pragma omp for
      for (int k = 1; k < nz-1; ++k) {
        for (int j = 1; j < ny-1; ++j) {
          for (int i = (k+j % 2); i < nx-1; i+=2) {
            int ix = i*dx + j*dy + k*dz;
            int _ix = 3*ix;
            half *_d = d + _ix;
            half *_b = b + _ix;
            half const *_sigma = sigma + ix;
            half *_phi = phi + ix;

            half p = half(1.0/6.0) * ( _sigma[0]
              + _phi[dx] + _phi[-dx] + _phi[dy] + _phi[-dy] + _phi[dx] + _phi[-dz]
              + _b[0] - _b[0-_dx] + _b[1] - _b[1-_dy] + _b[2] - _b[2-_dz]
              - _d[0] + _d[0-_dx] - _d[1] + _d[1-_dy] - _d[2] + _d[2-_dz] );

            _phi[0] = p > half(1.0) ? half(1.0) : (p < half(0.0) ? half(0.0) : p);
          }
        }
      } // update phi (even)

      // update phi (odd)
#pragma omp for
      for (int k = 1; k < nz-1; ++k) {
        for (int j = 1; j < ny-1; ++j) {
          for (int i = (k+j+1 % 2); i < nx-1; i+=2) {
            int ix = i*dx + j*dy + k*dz;
            int _ix = 3*ix;
            half *_d = d + _ix;
            half *_b = b + _ix;
            half const *_sigma = sigma + ix;
            half *_phi = phi + ix;

            half p = half(1.0/6.0) * ( _sigma[0]
              + _phi[dx] + _phi[-dx] + _phi[dy] + _phi[-dy] + _phi[dx] + _phi[-dz]
              + _b[0] - _b[0-_dx] + _b[1] - _b[1-_dy] + _b[2] - _b[2-_dz]
              - _d[0] + _d[0-_dx] - _d[1] + _d[1-_dy] - _d[2] + _d[2-_dz] );

            _phi[0] = p > half(1.0) ? half(1.0) : (p < half(0.0) ? half(0.0) : p);
          }
        }
      } // update phi (odd)

      // update b
#pragma omp for
      for (int k = 0; k < nz-1; ++k) {
        for (int j = 0; j < ny-1; ++j) {
          for (int i = 0; i < nx-1; ++i) {
            int ix = i*dx + j*dy + k*dz;
            int _ix = 3*ix;
            half *_d = d + _ix;
            half *_b = b + _ix;
            half *_phi = phi + ix;

            _b[0] += _phi[dx] - _phi[0] - _d[0];
            _b[1] += _phi[dy] - _phi[0] - _d[1];
            _b[2] += _phi[dz] - _phi[0] - _d[2];
          }
        }
      } // update b

    } // nits
  } // #pragma omp parallel
}
