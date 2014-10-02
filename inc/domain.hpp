#pragma once

namespace aminda {

  class Chunk {
  public:
    Chunk(int start, int stride, int length, int dx, int dy, int dz) : 
      _curr(start), _stride(stride), _end(start+length*stride),
      _dx(dx), _dy(dy), _dz(dz) {}
    void next() { _curr += _stride; }
    bool done() { _curr == _end; }
    int operator*() const { return _curr; }
    int px() const { return _curr + _dx; }
    int mx() const { return _curr - _dx; }
    int py() const { return _curr + _dy; }
    int my() const { return _curr - _dy; }
    int pz() const { return _curr + _dz; }
    int mz() const { return _curr - _dz; }
  private:
    int _curr;
    int _stride;
    int _end;
    int _dx;
    int _dy;
    int _dz;
  };

  class Domain {
  public:
    virtual int numAllButBackChunks() const = 0;
    virtual Chunk allButBackChunk(int i) const = 0;
    virtual int numInteriorEvenChunks() const = 0;
    virtual Chunk interiorEvenChunk(int i) const = 0;
    virtual int numInteriorOddChunks() const = 0;
    virtual Chunk interiorOddChunk(int i) const = 0;
  };

  class DenseDomain : public Domain {
  public:
    DenseDomain(int nx, int ny, int nz) : _nx(nx), _ny(ny), _nz(nz) {}

    int numAllButBackChunks() const { return (_ny-1)*(_nz-1); }
    int numInteriorEvenChunks() const { return (_ny-2)*(_nz-2); }
    int numInteriorOddChunks() const { return (_ny-2)*(_nz-2); }

    Chunk allButBackChunk(int ix) const {
      int y = ix%(_ny-1);
      int z = ix/(_ny-1);
      return Chunk(
    }

  private:
    const int _nx;
    const int _ny;
    const int _nz;
  };
}
