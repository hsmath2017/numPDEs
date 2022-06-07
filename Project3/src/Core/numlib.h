#ifndef NUMLIB_H
#define NUMLIB_H

#include "Vec.h"

/* 
  This file provides common functionalities in numerical computing. 
 
  1. common math functions.
  2. root-finding for non-linear equations.
  3. numerical quadrature.
*/

inline constexpr int sign(Real a) {
  if(a > 0) return 1;
  if(a < 0) return -1;
  return 0;
}

inline constexpr int factorial(int n) {
  int r = 1;
  for(; n>1; r*=n, --n);
  return r;
}

inline constexpr int binom(int n, int m) {
  return factorial(n) / factorial(n-m) / factorial(m);
}

template <class T, int p, int m = p%2>
struct ipow_wrapper;

template <class T, int p>
struct ipow_wrapper<T,p,0> {
  static constexpr auto calc(const T &a) {
    const auto x = ipow_wrapper<T,p/2>::calc(a);
    return x*x;
  }
};

template <class T, int p>
struct ipow_wrapper<T,p,1> {
  static constexpr auto calc(const T &a) {
    const auto x = ipow_wrapper<T,p/2>::calc(a);
    return x*x*a;
  }
};

template <class T>
struct ipow_wrapper<T,0,0> {
  static constexpr T calc(const T &) { return static_cast<T>(1); }
};

template <int p, class T> inline constexpr auto ipow(const T &a) {
  return ipow_wrapper<T,p>::calc(a);
}

//============================================================

template <class T_Func, class T_Der>
inline
void backtrack(const T_Func &f,
               const T_Der &df,
               Real &x,
               Real initialStep,           // the initial step size to go
               Real stepConstraint,       // the maximum step size to go
               int outerAttempts,         // number of iterations (if succesful) to take
               Real stepMultiplier = 1.1, // the factor by which the step size is multiplied if successful
               int innerAttempts = 10)    // number of iterations to take if the attempt is not successful
{
  Real fx = f(x);
  Real step = initialStep;
  while(outerAttempts-- > 0) {
    Real dfx = df(x);
    Real s = -sign(fx / dfx);
    bool successful = false;
    for(int i = 0; i < innerAttempts; ++i) {
      Real newx = x + s * step;
      Real newfx = f(newx);
      if(std::abs(newfx) < std::abs(fx)) {
        x = newx;
        fx = newfx;
        step = std::min(stepConstraint, step * stepMultiplier);
        successful = true;
        break;
      }
      step /= 2.0;
    }
    if(!successful)
      return;
  }
}

template <class T_Func, class T_Der>
inline
Real fzero(const T_Func &f, const T_Der &df, Real x0, int maxIter, Real tol)
{
  Real fx = f(x0);
  while(maxIter-- > 0 && std::abs(fx) > tol) {
    Real dfx = df(x0);
    x0 -= fx / dfx;
    fx = f(x0);
  }
  if(maxIter == 0 && std::abs(fx) > tol)
    dbgcout2 << "Newton iteration may not converge, f(x) = " << fx << std::endl;
  return x0;
}

//============================================================
    
template <int numK>
struct GaussLegendreConstant;

template <int numK, class T_Func>
inline Real quad(const T_Func &g, Real a, Real b) {
  using GLC = GaussLegendreConstant<numK>;
  Real R = 0;
  Real u = (a+b)/2, v = (b-a)/2;
  int k = 0;
  for(; k<(numK+1)/2; ++k)
    R += GLC::weights[k] * g(u + GLC::knots[k] * v);
  for(; k<numK; ++k)
    R += GLC::weights[numK-1-k] * g(u - GLC::knots[numK-1-k] * v);
  return v*R;
}

// Iterated integral on a multi-dimensional box

template <int Dim>
struct quad_wrapper
{
  template <int numK, class T_Func>
  static Real q(const T_Func &g, const Vec<Real, Dim> &lo, const Vec<Real, Dim> &hi) {
    auto lo_1 = reduce(lo);
    auto hi_1 = reduce(hi);
    auto inner = [&](Real t) {
      return quad_wrapper<Dim-1>::template q<numK>(
          [&](const Vec<Real, Dim-1> &x) { return g(enlarge(x, t)); },
          lo_1, hi_1);
    };
    return quad<numK>(inner, lo[Dim-1], hi[Dim-1]);
  }
};

template <>
struct quad_wrapper<1>
{
  template <int numK, class T_Func>
  static Real q(const T_Func &g, const Vec<Real, 1> &lo, const Vec<Real, 1> &hi) {
    return quad<numK>(g, lo[0], hi[0]);
  }
};

template <int numK, class T_Func, int Dim>
inline Real quad(const T_Func &g, const Vec<Real, Dim> &lo, const Vec<Real, Dim> &hi)
{
  return quad_wrapper<Dim>::template q<numK>(g, lo, hi);
}

// composite quadrature
template <class T_Func>
inline Real cquad(const T_Func &g, Real a, Real b, int subdiv) {
  Real R = 0.0;
  for(int k=0; k<subdiv; ++k)
    R += quad<3>(g, a + (b-a)*k/subdiv, a + (b-a)*(k+1)/subdiv);
  return R;
}

// adaptive Simpson's rule
template <class T_Func>
inline Real aquad(const T_Func &g, Real a, Real b, Real absTol = 1e-16) {
  const int minDepth = 1;
  const int maxDepth = 8;
  Real R[maxDepth+1][2]; // [][0] results of trapezoidal rule, [][1] results of Simpson's rule
  auto T = [&](int N) {
    Real z = 0.0;
    for(int i=0; i<N; ++i)
      z += g(a + (i+0.5)/N * (b-a));
    return z * (b-a) / N;
  };
  R[0][0] = (g(a) + g(b)) * (b-a) / 2;
  int k;
  Real delta;
  for(k=0; k<minDepth; ++k)
    R[k+1][0] = (R[k][0] + T(1<<k))/2;
  --k;
  R[k+1][1] = (4.0 * R[k+1][0] - R[k][0]) / 3.0;
  for(++k; k<maxDepth;) {
    R[k+1][0] = (R[k][0] + T(1<<k))/2;
    R[k+1][1] = (4.0 * R[k+1][0] - R[k][0]) / 3.0;
    delta = (R[k+1][1] - R[k][1]) / 15.0;
    ++k;
    if(std::abs(delta) < absTol)
      break;
  }
  if(std::abs(delta) >= absTol)
    dbgcout3 << "aquad() : est. error = " << delta << ", subdiv = " << (1<<k) << std::endl;
  return R[k][1] + delta;
}

//==================================================

template <>
struct GaussLegendreConstant<1>
{
  static constexpr Real knots[] = {0};
  static constexpr Real weights[] = {2.0};
};

template <>
struct GaussLegendreConstant<2>
{
  static constexpr Real knots[] = {-0.577350269189626};
  static constexpr Real weights[] = {1.0};
};

template <>
struct GaussLegendreConstant<3>
{
  static constexpr Real knots[] = {-0.774596669241483, 0};
  static constexpr Real weights[] = {5.0/9, 8.0/9};
};

template <>
struct GaussLegendreConstant<4>
{
  static constexpr Real knots[] = {-0.861136311594053, -0.339981043584856};
  static constexpr Real weights[] = {0.3478548451374538, 0.6521451548625461};
};

template <>
struct GaussLegendreConstant<5>
{
  static constexpr Real knots[] = {-0.906179845938664, -0.538469310105683, 0};
  static constexpr Real weights[] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889};
};

template <>
struct GaussLegendreConstant<6>
{
  static constexpr Real knots[] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969};
  static constexpr Real weights[] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910};
};

template <>
struct GaussLegendreConstant<7>
{
  static constexpr Real knots[] = {-0.9491079123427585, -0.7415311855993945,
                                   -0.4058451513773972, 0};
  static constexpr Real weights[] = {0.1294849661688697, 0.2797053914892766,
                                     0.3818300505051189, 0.4179591836734694};
};

//==================================================

namespace NS4_Ops {

template <class = void>
struct absolute {
  template <class T>
  constexpr auto operator() (const T &rhs) const { return std::abs(rhs); }
};

template <class = void>
struct max_of {
  template <class T>
  constexpr auto operator() (const T &lhs, const T &rhs) const { return (lhs < rhs) ? (rhs) : (lhs); }
};

}

#endif //NUMLIB_H
