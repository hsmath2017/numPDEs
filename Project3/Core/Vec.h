#ifndef VEC_H
#define VEC_H

#include <utility>
#include <initializer_list>
#include <cmath>
#include <cassert>
#include "Config.h"

// todo : support for constexpr

/// Define a point or a vector in the Euclidean space.
/**
 */
template <class T, int _Dim>
class Vec
{
protected:
  T coord[_Dim];

public:
  enum { Dim = _Dim };

  Vec(const T &t = T()) {
    for(int i=0; i<Dim; coord[i++] = t);
  }

  Vec(std::initializer_list<T> l) {
    auto j = l.begin();
    for(int d=0; d<Dim; ++d)
      coord[d] = *j++;
  }
  Vec(const Vec&) = default;
  Vec(Vec&&) = default;
  Vec& operator=(const Vec&) = default;
  Vec& operator=(Vec&&) = default;

  ///
  /**
    The conversion constructor.
   */
  template <class T2>
  explicit Vec(const Vec<T2,Dim>& rhs) {
    for(int d=0; d<Dim; d++)
      coord[d] = static_cast<T>(rhs[d]);
  }

  static Vec<T,Dim> unit(int D) {
    Vec<T,Dim> r;
    r[D] = static_cast<T>(1);
    return r;
  }

  // accessors
public:

  T &operator[](int _d) { return coord[_d]; }

  const T &operator[](int _d) const { return coord[_d];  }

  const T *data() const { return &coord[0]; }

public:

#define ELMWISE_BINARY_OP(OpName, Op)           \
  template <class T2>                           \
  auto OpName (const Vec<T2,Dim> &rhs) const {  \
    using Tx = decltype(coord[0] Op rhs[0]);    \
    Vec<Tx,Dim> res;                            \
    for(int i=0; i<Dim; i++)                    \
      res[i] = coord[i] Op rhs[i];              \
    return res;                                 \
  }

  ELMWISE_BINARY_OP(operator+, +)
  ELMWISE_BINARY_OP(operator-, -)
  ELMWISE_BINARY_OP(operator*, *)
  ELMWISE_BINARY_OP(operator/, /)
#undef ELMWISE_BINARY_OP

#define RIGHT_BROADCAST(OpName, Op)        \
  template <class T2>                      \
  auto OpName (const T2 &rhs) const {      \
    using Tx = decltype(coord[0] Op rhs);  \
    Vec<Tx,Dim> res;                       \
    for(int d=0; d<Dim; ++d)               \
      res[d] = coord[d] Op rhs;            \
    return res;                            \
 }

  RIGHT_BROADCAST(operator+, +)
  RIGHT_BROADCAST(operator-, -)
  RIGHT_BROADCAST(operator*, *)
  RIGHT_BROADCAST(operator/, /)
#undef RIGHT_BROADCAST

  ///
  /**
     Formatting.
   */
  friend std::ostream &operator<<(std::ostream &os, const Vec<T,Dim> &p) {
    os << "(" << p[0];
    for(int d=1; d<Dim; d++)
      os << "," << p[d];
    os << ")";
    return os;
  }
};

//==================================================

/// 2D cross product gives a scalar
template <class T>
inline T cross(const Vec<T,2> &va, const Vec<T,2> &vb)
{
  return va[0]*vb[1] - va[1]*vb[0];
}

/// 3D cross product gives a vector
template <class T>
inline Vec<T,3> cross(const Vec<T,3> &va, const Vec<T,3> &vb)
{
  return Vec<T,3> {
      va[1]*vb[2] - va[2]*vb[1],
      va[2]*vb[0] - va[0]*vb[2],
      va[0]*vb[1] - va[1]*vb[0]
  };
}

template <class T>
inline Vec<T,2> clockwise(const Vec<T,2> &v) {
  return Vec<T,2> {v[1], -v[0]};
}

//==================================================
// operations on multi-dim vectors

template <class T, int Dim>
inline Vec<T,Dim> operator-(const Vec<T,Dim> &lhs) {
  Vec<T,Dim> res;
  for(int d=0; d<Dim; d++)
    res[d] = -lhs[d];
  return res;
}

template <int Dim>
inline Vec<Real,Dim> floor(const Vec<Real,Dim> &rhs) {
  Vec<Real,Dim> res;
  for(int d=0; d<Dim; res[d] = std::floor(rhs[d]), ++d);
  return res;
}

template <int Dim>
inline Vec<Real,Dim> abs(const Vec<Real,Dim> &rhs) {
  Vec<Real,Dim> res;
  for(int d=0; d<Dim; res[d] = std::abs(rhs[d]), ++d);
  return res;
}

template <class T, int Dim>
inline Vec<int,Dim> sign(const Vec<T,Dim>& a) {
  Vec<int,Dim> res;
  for(int d=0; d<Dim; d++)
    res[d] = ((a[d] > 0) ? (1) : ((a[d] < 0) ? (-1) : (0)));
  return res;
}

inline Real norm(Real x, int nt)
{
  return std::abs(x);
}

template <class T, int Dim>
inline T norm(const Vec<T,Dim> &lhs, int nt = 2)
{
  if(nt == 2) {
    return sqrt(dot(lhs, lhs));
  } else if(nt == 1) {
    T a = 0;
    for(int i=0; i<Dim; a+=std::abs(lhs[i]), ++i);
    return a;
  } else if(nt == 0) {
    T a = std::abs(lhs[0]);
    for(int i=1; i<Dim; i++)
      a = std::max(std::abs(lhs[i]), a);
    return a;
  }
  return 0; // never reached
}

template <int Dim>
inline Vec<Real,Dim> normalize(const Vec<Real,Dim> &lhs) {
  Real l = norm(lhs, 2);
  return lhs / l;
}

template <class T, int Dim>
inline T dot(const Vec<T,Dim> &lhs, const Vec<T,Dim> &rhs) {
  T r = 0;
  for(int d=0; d<Dim; d++)
    r += lhs[d] * rhs[d];
  return r;
}

template <class T, int Dim>
inline T sum(const Vec<T,Dim> &v)
{
  auto a = v[0];
  for(int d=1; d<Dim; a+=v[d], ++d);
  return a;
}

template <class T, int Dim>
inline T prod(const Vec<T,Dim> &v)
{
  auto a = v[0];
  for(int d=1; d<Dim; a*=v[d], ++d);
  return a;
}

template <class T, int Dim, class Cmp = std::less<T>>
inline T min_of(const Vec<T,Dim> &v, Cmp cmp = Cmp())
{
  auto a = v[0];
  for(int d=1; d<Dim; d++)
    if(cmp(v[d], a))
      a = v[d];
  return a;
}

template <class T, int Dim>
inline T max_of(const Vec<T,Dim>& v)
{
  return min_of(v, std::greater<T>());
}

// element-wise min
template <class T, int Dim>
inline Vec<T,Dim> min(const Vec<T,Dim>& lhs, const Vec<T,Dim>& rhs)
{
  Vec<T,Dim> res;
  for(int d=0; d<Dim; d++)
    res[d] = std::min(lhs[d], rhs[d]);
  return res;
}

// element-wise max
template <class T, int Dim>
inline
Vec<T,Dim> max(const Vec<T,Dim>& lhs, const Vec<T,Dim>& rhs)
{
  Vec<T,Dim> res;
  for(int d=0; d<Dim; d++)
    res[d] = std::max(lhs[d], rhs[d]);
  return res;
}

//==================================================

template <class T, int Dim>
inline Vec<T,Dim-1> reduce(const Vec<T,Dim> &lhs, int D = Dim-1) {
  Vec<T,Dim-1> rhs;
  int i;
  for(i=0; i<D; ++i)
    rhs[i] = lhs[i];
  for(; i<Dim-1; ++i)
    rhs[i] = lhs[i+1];
  return rhs;
}

template <class T, int Dim>
inline Vec<T,Dim+1> enlarge(const Vec<T,Dim> &lhs, const T &a, int k = Dim) {
  Vec<T,Dim+1> res;
  for(int d=0; d<k; res[d] = lhs[d], d++);
  res[k] = a;
  for(int d=k; d<Dim; res[d+1] = lhs[d], d++);
  return res;
}

#endif // VEC_H
