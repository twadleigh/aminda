#pragma once

#include <double_conv.hpp>
#include <vec.hpp>
#include <domain.hpp>

namespace aminda {

  namespace _internal {

    template<typename O, typename U>
    inline void applyWithGradient(const O& op, const Domain& dom, const U* u) {
      const int NC = dom.numAllButBackChunks();
#pragma omp parallel for
      for(int k = 0; k < NC; ++k) {
	for(Chunk c = dom.allButBackChunk(k); !c.done(); c.next()) {
	  op(c0,vec<3,double>(toDouble<U>(u[c0])-toDouble<U>(u[c.mx()]),
			      toDouble<U>(u[c0])-toDouble<U>(u[c.my()]),
			      toDouble<U>(u[c0])-toDouble<U>(u[c.mz()])));
	}
      }
    }
  }
  
  template<typename B, typename U, typename D>
  void updateB(const Domain& dom, vec<3,B>* b, const U* u, const vec<3,D>* d) {
    struct O {
      vec<3,B>* _b;
      const vec<3,D>* _d;
      O(vec<3,B>* bb, const vec<3,D>* dd) : _b(bb), _d(dd) {}
      void operator()(int c0, vec<3,double> gradU) const {
	_b[c0][0] += fromDouble<B>(gradU[0]-toDouble<D>(_d[c0][0]));
	_b[c0][1] += fromDouble<B>(gradU[1]-toDouble<D>(_d[c0][1]));
	_b[c0][2] += fromDouble<B>(gradU[2]-toDouble<D>(_d[c0][2]));
      }
    } op(b,d);
    _internal::applyWithGradient<O,U>(op,dom,u);
  }

  namespace _internal {
    template<typename R, typename D, typename U, typename B>
    void updateD(const R& r, const Domain& dom, vec<3,D>* d, const U* u, 
		 const vec<3,B>* b) {
      struct O {
	vec<3,D>* _d;
	const vec<3,B>* _b;
	const R& _r;
	O(vec<3,D>* dd, const vec<3,B>* bb, const R& rr) : 
	  _d(dd), _b(bb), _r(rr) {}
	void operator()(int c0, vec<3,double> gradU) const {
	  vec<3,double> newD = 
	    vec<3,double>(gradU[0]+toDouble<B>(_b[c0][0]),
			  gradU[1]+toDouble<B>(_b[c0][1]),
			  gradU[2]+toDouble<B>(_b[c0][2]));
	  double lenD = sqrt(newD[0]*newD[0]+newD[1]*newD[1]+newD[2]*newD[2]);
	  double thresh = _r(c0);
	  double scale = lenD > thresh ? 1.0-(thresh/lenD);
	  _d[c0][0] = fromDouble<D>(scale*newD[0]);
	  _d[c0][1] = fromDouble<D>(scale*newD[1]);
	  _d[c0][2] = fromDouble<D>(scale*newD[2]);
	}
      } op(d,b,r);
      _internal::applyWithGradient<O,U>(op,dom,u);
    }
  }

  template<typename D, typename U, typename B, typename R>
  void updateD(const Domain& dom, vec<3,D>* d, const U* u, const vec<3,B>* b, 
	       const R* r, double l) {
    struct RR {
      const R* _r;
      const double _l;
      RR(const R* rr, double ll) : _r(rr), _l(ll) {}
      double operator()(int c0) const { return _r[c0]*_l; }
    } rr(r,l);
    _internal::updateD<RR,D,U,B>(rr,dom,d,u,b);
  }

  template<typename D, typename U, typename B>
  void updateDUnitR(const Domain& dom, vec<3,D>* d, const U* u, const vec<3,B>* b, 
		    double l) {
    struct RR {
      const double _l;
      RR(double ll) : _l(ll) {}
      double operator()(int) const { return _l; }
    } rr(l);
    _internal::updateD<RR,D,U,B>(rr,dom,d,u,b);
  }

  namespace _internal {
    template<typename S, typename U, typename D, typename B>
    inline void updateUChunk(const S& s, Chunk& c, U* u, const vec<3,D>* d, 
		      const vec<3,B>* b) {
      const double oneSixth = 1.0/6.0;
      for(; !c.done(); c.next()) {
	const int c0 = *c;
	const int cx = c.px();
	const int cy = c.py();
	const int cz = c.pz();
	double newU = oneSixth*
	  (toDouble<U>(u[cx])+toDouble<U>(u[c.mx()])+
	   toDouble<U>(u[cy])+toDouble<U>(u[c.my()])+
	   toDouble<U>(u[cz])+toDouble<U>(u[c.mz()])+
	   toDouble<B>(b[cx][0])-toDouble<B>(b[c0][0])-
	   toDouble<D>(d[cx][0])+toDouble<D>(d[c0][0])+
	   toDouble<B>(b[cy][1])-toDouble<B>(b[c0][1])-
	   toDouble<D>(d[cy][1])+toDouble<D>(d[c0][1])+
	   toDouble<B>(b[cz][2])-toDouble<B>(b[c0][2])-
	   toDouble<D>(d[cz][2])+toDouble<D>(d[c0][2])-
	   s(c0));
	u[c0] = fromDouble<U>(newU > 0.0 ? (newU < 1.0 ? newU : 1.0) : 0.0);
      }
    }

    template<typename S, typename U, typename D, typename B>
    inline void updateU(const S& s, const Domain& dom, U* u, const vec<3,D>* d,
		 const vec<3,B>* b) {
      const int NE = dom.numInteriorEvenChunks();
#pragma omp parallel for
      for(int k = 0; k < NE; ++k) {
	updateUChunk<S,U,D,B>(s,dom.interiorEvenChunk(k),u,d,b);
      }
      const int NO = dom.numInteriorOddChunks();
#pragma omp parallel for
      for(int k = 0; k < NO; ++k) {
	updateUChunk<S,U,D,B>(s,dom.interiorOddChunk(k),u,d,b);
      }
  }

  template<typename U, typename D, typename B, typename S>
  void updateU(const Domain& dom, U* u, const vec<3,D>* d, const vec<3,B>* b,
	       const S* s, double l) {
    struct SS {
      const S* _s;
      const double _l;
      SS(const S* ss, double ll) : _s(ss), _l(ll) {}
      double operator()(int c0) const {	return _l*toDouble<S>(s[c0]); }
    } sigma(s,l);
    _internal::updateU<SS,U,D,B>(sigma,dom,u,d,b);
  }

  template<typename U, typename D, typename B>
  void updateUUnitS(const Domain& dom, U* u, const vec<3,D>* d, const vec<3,B>* b,
		    double l) {
    struct SS {
      const double _l;
      SS(double ll) : _l(ll) {}
      double operator()(int) const { return _l*toDouble<S>(s[c0]); }
    } sigma(l);
    _internal::updateU<SS,U,D,B>(sigma,dom,u,d,b);
  }

  template<typename U, typename D, typename B>
  void updateUZeroS(const Domain& dom, U* u, const vec<3,D>* d, const vec<3,B>* b)
    struct SS {
      double operator()(int) const { return 0.0; }
    } sigma;
    _internal::updateU<SS,U,D,B>(sigma,dom,u,d,b);
  }

}
