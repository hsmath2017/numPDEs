#ifndef INTERVAL_H
#define INTERVAL_H

#include "Vec.h"

///
/**
   Represent a box in \mathbb{R}^D.
 */
template <int Dim>
class Interval
{
public:
  using rVec = Vec<Real,Dim>;

  Interval() = default;

  Interval(const rVec &_lo, const rVec &_hi) : corners{_lo, _hi} { }

  const rVec &lo() const { return corners[0]; }
  const rVec &hi() const { return corners[1]; }

  // arithmetics
public:
  Interval<Dim> operator+(const rVec &x) const {
    return Interval<Dim>(corners[0]+x, corners[1]+x);
  }
  Interval<Dim> operator-(const rVec &x) const {
    return Interval<Dim>(corners[0]-x, corners[1]-x);
  }
  friend Interval<Dim> operator*(Real s, const Interval<Dim> &itv) {
    return Interval<Dim>(itv.lo()*s, s*itv.hi()*s);
  }
  friend Interval<Dim> operator/(const Interval<Dim> &itv, Real s) {
    return Interval<Dim>(itv.lo()/s, itv.hi()/s);
  }

public:
  ///
  /**
     Intersection.
   */
  Interval<Dim> operator& (const Interval<Dim> &rhs) const {
    return Interval<Dim>(max(lo(), rhs.lo()),
                         min(hi(), rhs.hi()));
  }

  bool contain(const rVec& t, Real tol = 0) const {
    return (min_of(hi()-t) >= -tol) && (min_of(t-lo()) >= -tol);
  }

  bool contain(const Interval<Dim>& rhs, Real tol = 0) const {
    return contain(rhs.lo(), tol) && contain(rhs.hi(), tol);
  }

  bool empty(Real tol = 0) const {
    return min_of(hi()-lo()) <= tol;
  }

  Real volume() const { return prod(hi()-lo()); }

protected:
  rVec corners[2];
};

template <int Dim>
inline Interval<1> getComp(const Interval<Dim> &itv, int comp) {
  return Interval<1>(itv.lo()[comp], itv.hi()[comp]);
}

#endif // INTERVAL_H
