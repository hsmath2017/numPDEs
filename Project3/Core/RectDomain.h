#ifndef RECTGRID_H
#define RECTGRID_H

#include <array>
#include "Box.h"

enum { FaceCentered = 0/*From 0 to Dim-1*/, CellCentered = -1, NodeCentered = -2 };

template <int Dim>
class RectDomain : public Box<Dim>
{
public:
  using rVec = Vec<Real,Dim>;
  using iVec = Vec<int,Dim>;
  using BaseClass = Box<Dim>;

  RectDomain() = default;

  RectDomain(const BaseClass &aBox, const rVec& aDx, int aCentering, int aNumGhost)
      : BaseClass(aBox), dx(aDx), centering(aCentering), nGhost(aNumGhost)
  {
  }

  // accessor
public:
  using BaseClass::lo;
  using BaseClass::hi;
  using BaseClass::size;
  using BaseClass::volume;
  //
  const rVec &spacing() const { return dx; }
  const rVec getDelta() const {
    rVec delta = 0.5;
    if(centering == CellCentered) {
      return delta;
    } else if(centering == NodeCentered) {
      return rVec(0.0);
    } else {
      delta[centering] = 0;
      return delta;
    }
  }
  int getCentering() const { return centering; }
  int getNumGhost() const { return nGhost; }
  //
  BaseClass getGhostedBox() const { return BaseClass::inflate(nGhost); }

  // deformation
public:
  RectDomain refine() const {
    iVec newlo = lo() * 2;
    iVec newhi;
    if(centering == NodeCentered) {
      newhi = hi() * 2;
    } else if(centering == CellCentered) {
      newhi = hi() * 2 + 1;
    } else if(centering >= 0) {
      newhi = hi() * 2 + 1 - iVec::unit(centering);
    }
    return RectDomain(BaseClass(newlo, newhi), dx/2, centering, nGhost);
  }
  RectDomain coarsen() const {
    iVec newlo = lo() / 2;
    iVec newhi = hi() / 2;
    return RectDomain(BaseClass(newlo, newhi), dx*2, centering, nGhost);
  }
  RectDomain stagger(int type) const {
    if(type == centering)
      return *this;
    Box<Dim> ccBox;
    if(centering == CellCentered) {
      ccBox = *this;
    } else if(centering == NodeCentered) {
      ccBox = Box<Dim>(lo(), hi()-1);
    } else {
      ccBox = Box<Dim>(lo(), hi() - iVec::unit(centering));
    }
    if(type == CellCentered) {
      return RectDomain(ccBox, dx, type, nGhost);
    } else if(type == NodeCentered) {
      return RectDomain({ccBox.lo(), ccBox.hi()+1}, dx, type, nGhost);
    } else {
      return RectDomain({ccBox.lo(), ccBox.hi() + iVec::unit(type)}, dx, type, nGhost);
    }
  }

protected:
  rVec dx;
  int centering;
  int nGhost;
};

#endif // RECTGRID_H
