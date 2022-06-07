#ifndef CURVE_H
#define CURVE_H

#include <vector>
#include "Polynomial.h"


template <int Dim, int Order>
class Curve
{
 public:
  // fitCurve's boundary type
  enum BCType {
    notAknot = 0,
    periodic,
    complete,
    second,
    nature,
    nBC_type
  };
  using rVec = Vec<Real,Dim>;
  using T_Polynomial = Polynomial<Order,rVec>;
  template <class T> using vector = std::vector<T>;

protected:
  std::vector<Real> knots;
  std::vector<Polynomial<Order,rVec>> polys;

public:
  // constructors
  Curve() = default;

  Curve(const Polynomial<Order,rVec>& pn, Real len) {
    polys.push_back(pn);
    knots.push_back(0);
    knots.push_back(len);
  }

  Curve(const Curve&) = default;
  auto operator=(const Curve&) -> Curve& = default;
  Curve(Curve&&) = default;
  auto operator=(Curve&&) -> Curve& = default;

  ~Curve() = default;

  // accessors
public:
  // helpers
  int locatePiece(Real t) const;

  rVec operator() (Real t) const {
    int p = locatePiece(t);
    return (polys[p])(t - knots[p]);
  }

  Interval<1> domain() const {
    if(empty())
      return Interval<1>(0.0, 0.0);
    return Interval<1>(knots.front(), knots.back());
  }

  bool empty() const {
    return knots.size() <= 1;
  }

  rVec startpoint() const {
    return polys.front()[0];
  }

  rVec endpoint() const {
    auto np = polys.size();
    Real t = knots[np] - knots[np-1];
    return polys[np-1](t);
  }

  rVec midpoint() const {
    return (*this)(0.5*(knots.front() + knots.back()));
  }

  bool isClosed(Real tol) const {
    return norm(endpoint() - startpoint(), 2) <= tol;
  }

  const std::vector<Real> &getKnots() const { return knots; }

  const std::vector<T_Polynomial> &getPolys() const { return polys; }

  // modifiers
public:

  void concat(const T_Polynomial &p, Real plen);

  void concat(const Curve<Dim,Order> &pp);

  Curve<Dim,Order> reverse() const;

  Curve<Dim,Order> offset(const rVec &ofs) const {
    Curve<Dim,Order> res;
    res.knots = this->knots;
    for(const auto& p : this->polys)
      res.polys.push_back(p + ofs);
    return res;
  }

  Curve<Dim,Order> rotate(Real theta) const {
    const Real c = cos(theta), s = sin(theta);
    Curve<Dim, Order> res;
    res.knots = this->knots;
    for(const auto &p : this->polys) {
      T_Polynomial q;
      for(int k = 0; k < Order; ++k) {
        q[k][0] = p[k][0] * c + p[k][1] * (-s);
        q[k][1] = p[k][0] * s + p[k][1] * c;
      }
      res.polys.push_back(q);
    }
    return res;
  }

  // advanved operations
public:
  Curve<Dim,Order> extract(Real from, Real to, Real tol) const;

  void split(const vector<Real> &brks, vector<Curve<Dim,Order>> &out, Real tol) const;
  ///
  /**
     Return a copy of this curve,
     so that each piece is monotonic in both x and y direction.
   */
  Curve<Dim,Order> makeMonotonic(Real tol) const;
  ///
  /**
     Report the number of proper intersections with the line x_d=c.
     Require the curve to be piecewise monotonic, see makeMonotonic().
   */
  int countProperInts(Real c, int d, Real tol) const;

  // friend helpers
public:
  template <int Ord>
  friend Curve<2,Ord> createRect(const Vec<Real,2>& lo, const Vec<Real,2>& hi);
  template <int Ord>
  friend Curve<2,Ord> createLineSegment(const Vec<Real,2>& p0, const Vec<Real,2>& p1);
  template <int Dm, int Ord>
  friend Curve<Dm,Ord-1> der(const Curve<Dm,Ord> &c);
  // template <int Ord>
  // friend Curve<2,Ord> fitCurve(const std::vector<Vec<Real,2>> &knots,
  //                              bool periodic);
  template <int Ord>
  friend Curve<2, Ord> fitCurve(const std::vector<Vec<Real, 2>>& knots,
                                typename Curve<2, Ord>::BCType type,
                                const Vec<Real, 2>& start,
                                const Vec<Real, 2>& end);
  // read/write operations
public:
  static Curve<Dim,Order> load(std::istream &is);
  void dump(std::ostream &os) const;
};

//================================================

template <int Dim, int Order>
Curve<Dim,Order-1> der(const Curve<Dim,Order> &c);

template <int Dim, int Order>
Interval<Dim> boundingBox(const Curve<Dim,Order> &c);

template <int Dim, int Order>
Interval<Dim> boundingBox(const std::vector<Curve<Dim,Order>> &vc);

template <int Order>
Real area(const Curve<2,Order> &gon);

template <int Order>
Real arclength(const Curve<2,Order> &c);

template <int Ord>
Curve<2,Ord> createRect(const Vec<Real,2>& lo, const Vec<Real,2>& hi);

template <int Ord>
Curve<2,Ord> createLineSegment(const Vec<Real,2>& p0, const Vec<Real,2>& p1);

template <int Order>
Curve<2,Order> fitCurve(const std::vector<Vec<Real,2>> &knots, typename Curve<2, Order>::BCType
                        type = Curve<2, Order>::notAknot,const Vec<Real,2>& start =
                        Vec<Real,2>(), const Vec<Real,2>& end = Vec<Real,2>());

#endif // CURVE_H
