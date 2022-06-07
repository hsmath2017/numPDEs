#ifndef BOX_H
#define BOX_H

#include "Vec.h"

template <int Dim>
class Box
{
public:
  using iVec = Vec <int,Dim>;

  Box() = default;

  Box(const iVec& _lo, const iVec& _hi) : corners{_lo, _hi} { }

  const iVec& lo() const { return corners[0]; }
  const iVec& hi() const { return corners[1]; }

  iVec size() const { return hi() - lo() + 1; }

  // deformation
public:
  Box<Dim> operator+(const iVec &delta) const {
    return Box<Dim>(lo()+delta, hi()+delta);
  }
  Box<Dim> operator-(const iVec &delta) const {
    return Box<Dim>(lo()-delta, hi()-delta);
  }

  // positive delta for inflation, negative delta for shrinkage
  Box<Dim> inflate(const iVec &delta) const {
    return Box<Dim>(lo()-delta, hi()+delta);
  }

  Box<Dim> operator & (const Box& rhs) const {
    auto l = max(lo(), rhs.lo());
    auto h = min(hi(), rhs.hi());
    return Box<Dim>(l,h);
  }

//  Box<Dim> refine() const { return Box<Dim>(corners[0]*2, corners[1]*2+1); }
//  Box<Dim> coarsen() const { return Box<Dim>(corners[0]/2, corners[1]/2); }

public:
  bool empty() const {
    return min_of(hi()-lo()) < 0;
  }

  bool contain(const iVec& pos) const {
    return (min_of(hi()-pos) >= 0) && (min_of(pos-lo()) >= 0);
  }

  bool contain(const Box<Dim>& rhs) const {
    return contain(rhs.lo()) && contain(rhs.hi());
  }

  int volume() const { return prod(size()); }

protected:
  iVec corners[2];
};

template <int Dim>
inline Box<1> getComp(const Box<Dim> &bx, int comp) {
  return Box<1>(bx.lo()[comp], bx.hi()[comp]);
}

//==================================================

template <int Dim>
inline Box<Dim+1> enlarge(const Box<Dim> &bx, const Box<1> &r, int k = Dim) {
  return Box<Dim+1>(enlarge(bx.lo(), r.lo()[0], k), enlarge(bx.hi(), r.hi()[0], k));
}

template <int Dim>
inline Box<Dim-1> reduce(const Box<Dim> &bx, int D = Dim-1)
{
  return Box<Dim-1>(reduce(bx.lo(), D), reduce(bx.hi(), D));
}

//==================================================

#define loop_box_1(bx, i0) \
for(int i0=bx.lo()[0]; i0<=bx.hi()[0]; ++i0)

#define loop_sz_1(sz, i0) \
for(int i0=0; i0<sz[0]; ++i0)

#define loop_box_2(bx, i0, i1) \
for(int i1=bx.lo()[1]; i1<=bx.hi()[1]; ++i1) \
for(int i0=bx.lo()[0]; i0<=bx.hi()[0]; ++i0)

#define loop_sz_2(sz, i0, i1) \
for(int i1=0; i1<sz[1]; ++i1) \
for(int i0=0; i0<sz[0]; ++i0)

#define loop_box_3(bx, i0, i1, i2) \
for(int i2=bx.lo()[2]; i2<=bx.hi()[2]; ++i2) \
for(int i1=bx.lo()[1]; i1<=bx.hi()[1]; ++i1) \
for(int i0=bx.lo()[0]; i0<=bx.hi()[0]; ++i0)

#define loop_sz_3(sz, i0, i1, i2) \
for(int i2=0; i2<sz[2]; ++i2) \
for(int i1=0; i1<sz[1]; ++i1) \
for(int i0=0; i0<sz[0]; ++i0)

#define loop_box_4(bx, i0, i1, i2, i3) \
for(int i3=bx.lo()[3]; i3<=bx.hi()[3]; ++i3) \
for(int i2=bx.lo()[2]; i2<=bx.hi()[2]; ++i2) \
for(int i1=bx.lo()[1]; i1<=bx.hi()[1]; ++i1) \
for(int i0=bx.lo()[0]; i0<=bx.hi()[0]; ++i0)

#define loop_sz_4(sz, i0, i1, i2, i3) \
for(int i3=0; i3<sz[3]; ++i3) \
for(int i2=0; i2<sz[2]; ++i2) \
for(int i1=0; i1<sz[1]; ++i1) \
for(int i0=0; i0<sz[0]; ++i0)

#endif //BOX_H
