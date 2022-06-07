#ifndef VECCOMPARE_H
#define VECCOMPARE_H

#include "Vec.h"

/// Give lexicological order on vectors
/**
 */
template <class T, int Dim>
class VecCompare;

///
/**
   Y coordinate as the major keyword (descending), 
   and X coordinate as the minor keyword (ascending).
 */
template <>
class VecCompare<Real,2>
{
public:
  enum { Dim=2 };
  typedef Vec<Real,Dim> rVec;
  
  ///
  /**
     The constructor accepts an absolute tolerance.
   */
  VecCompare(Real _tol = 0) : tol(_tol) { }
  
  ///
  /**
     Return the order between p1 and p2. 

     -1 for p1<p2, 0 for p1==p2, 1 for p1>p2. 

     Two points are considered equal if their distance 
     is about or less than the scale of tol.
  */
  int compare(const rVec &p1, const rVec &p2) const {
    if(std::abs(p1[1]-p2[1]) <= tol) {
      if(std::abs(p1[0]-p2[0]) <= tol) return 0;
      return (p1[0]<p2[0])?(-1):(1);
    }
    return (p1[1]>p2[1])?(-1):(1);
  }
  
  ///
  /**
     Return true if p1<p2. 
     It is a binary operator subject to strict weak ordering.
   */
  bool operator() (const rVec &p1, const rVec &p2) const
  {
    return (compare(p1,p2) == -1);
  }

protected:
  Real tol;
};

//==================================================
template <int Dim>
class VecCompare<int,Dim>
{
public:
  using iVec = Vec<int,Dim>;
  VecCompare() = default;

  // -1 for lhs<rhs, 0 for lhs==rhs, 1 for lhs>rhs.
  int compare(const iVec &lhs, const iVec &rhs) const {
    for(int d=Dim-1; d>=0; --d) {
      if(lhs[d] < rhs[d])
        return -1;
      if(lhs[d] > rhs[d])
        return 1;
    }
    return 0;
  }

  bool operator() (const iVec &lhs, const iVec &rhs) const {
    return compare(lhs, rhs) == -1;
  }
};

#endif // VECCOMPARE_H
