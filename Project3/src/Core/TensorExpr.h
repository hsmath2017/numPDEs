#ifndef TENSOREXPR_H
#define TENSOREXPR_H

#include <type_traits>
#include <limits>
#include <functional>
#include "numlib.h"
#include "Tensor.h"

// complete the definition of assignment from expression in Tensor<>
template <class T, int Dim>
template <class TExpr>
inline
Tensor<T,Dim>& Tensor<T,Dim>::operator=(const TensorExpr<TExpr> &expr)
{
  assert(norm(size() - expr.box().size(), 0) == 0);
  auto sz = expr.box().size();
  if(Dim == 1) {
#pragma omp parallel for default(shared) schedule(static)
    loop_sz_1(sz, i)
      at(i) = expr.at(i);
  } else if(Dim == 2) {
#pragma omp parallel for default(shared) schedule(static)
    loop_sz_2(sz, i, j)
      at(i, j) = expr.at(i, j);
  } else if(Dim == 3) {
#pragma omp parallel for default(shared) schedule(static)
    loop_sz_3(sz, i, j, k)
      at(i, j, k) = expr.at(i, j, k);
  } else if(Dim == 4) {
#pragma omp parallel for default(shared) schedule(static)
    loop_sz_4(sz, i, j, k, l)
      at(i, j, k, l) = expr.at(i, j, k, l);
  }
  return *this;
}

template <class T, int Dim>
template <class TExpr>
inline
Tensor<T,Dim>::Tensor(const TensorExpr<TExpr> &expr) : Tensor(expr.box())
{
  *this = expr;
}

//================================================

// if an expression does not carry box information,
// just return null_box_tag
struct null_box_tag { };

template <class T>
struct TensorLiteral : TensorExpr<TensorLiteral<T>>
{
  T a;
  TensorLiteral(const T &_a) : a(_a) { }

  template <class... Ts>
  T    at(int, Ts...) const { return a; }
  auto box()          const { return null_box_tag(); }
};

template <class TExpr, class Op>
struct TensorUnaryOp : public TensorExpr<TensorUnaryOp<TExpr,Op>>
{
  TExpr e1;
  TensorUnaryOp(const TExpr &e) : e1(e) { }

  template <class... Ts>
  auto        at(int i, Ts... args) const { return Op()(e1.at(i, args...)); }
  const auto  box()                 const { return e1.box(); }
};

template <class TExpr1, class TExpr2, class Op>
struct TensorBinaryOp : public TensorExpr<TensorBinaryOp<TExpr1,TExpr2,Op>>
{
  TExpr1 e1;
  TExpr2 e2;
  TensorBinaryOp(const TExpr1 &i1, const TExpr2 &i2) : e1(i1), e2(i2) { }

  template <class... Ts>
  auto        at(int i, Ts... args) const { return Op()(e1.at(i, args...), e2.at(i, args...)); }
  const auto  box()                 const { return box_by_tag(e1.box()); }

  template <int Dim> const auto& box_by_tag(const Box<Dim> &bx) const { return bx; }
  const auto&                    box_by_tag(null_box_tag)       const { return e2.box(); }
};

//==================================================
// unary
template <class TExpr>
inline
auto abs(const TensorExpr<TExpr> &rhs)
{
  return TensorUnaryOp<TExpr,NS4_Ops::absolute<>>(static_cast<const TExpr &>(rhs));
}

template <class TExpr>
inline
auto operator-(const TensorExpr<TExpr> &rhs)
{
  return TensorUnaryOp<TExpr,std::negate<>>(static_cast<const TExpr &>(rhs));
}

//==================================================
// both operands are expressions

template <class TExpr1, class TExpr2>
inline
auto operator+(const TensorExpr<TExpr1> &lhs, const TensorExpr<TExpr2> &rhs)
{
  return TensorBinaryOp<TExpr1,TExpr2,std::plus<>>
      (static_cast<const TExpr1 &>(lhs), static_cast<const TExpr2 &>(rhs)); // downcast to derived class
}

template <class TExpr1>
inline
auto operator+(const TensorExpr<TExpr1> &lhs, Real rhs)
{
  return TensorBinaryOp<TExpr1,TensorLiteral<Real>,std::plus<>>
      (static_cast<const TExpr1 &>(lhs), rhs); // downcast to derived class
}

template <class TExpr1, class TExpr2>
inline
auto operator-(const TensorExpr<TExpr1> &lhs, const TensorExpr<TExpr2> &rhs)
{
  return TensorBinaryOp<TExpr1,TExpr2,std::minus<>>
      (static_cast<const TExpr1 &>(lhs), static_cast<const TExpr2 &>(rhs)); // downcast to derived class
}

template <class TExpr1>
inline
auto operator-(const TensorExpr<TExpr1> &lhs, Real rhs)
{
  return TensorBinaryOp<TExpr1,TensorLiteral<Real>,std::minus<>>
      (static_cast<const TExpr1 &>(lhs), rhs); // downcast to derived class
}

template <class TExpr1, class TExpr2>
inline
auto operator*(const TensorExpr<TExpr1> &lhs, const TensorExpr<TExpr2> &rhs)
{
  return TensorBinaryOp<TExpr1,TExpr2,std::multiplies<>>
      (static_cast<const TExpr1 &>(lhs), static_cast<const TExpr2 &>(rhs)); // downcast to derived class
}

template <class TExpr1>
inline
auto operator*(const TensorExpr<TExpr1> &lhs, Real rhs)
{
  return TensorBinaryOp<TExpr1,TensorLiteral<Real>,std::multiplies<>>
      (static_cast<const TExpr1 &>(lhs), rhs);
}

//================================================
// other useful ops

template <class TExpr, class RedcOp, class TVal>
inline auto reduce(const TensorExpr<TExpr> &expr, const RedcOp &rop, const TVal &initVal)
{
  using T = std::decay_t<decltype(expr.at(0))>;
  T r = static_cast<T>(initVal);
  auto sz = expr.box().size();
  constexpr int Dim = decltype(sz)::Dim;
#pragma omp parallel default(shared)
  {
    T lr = static_cast<T>(initVal);
    if (Dim == 1) {
#pragma omp for schedule(static)
      loop_sz_1(sz, i)
        lr = rop(lr, expr.at(i));
    } else if (Dim == 2) {
#pragma omp for schedule(static)
      loop_sz_2(sz, i, j)
        lr = rop(lr, expr.at(i, j));
    } else if (Dim == 3) {
#pragma omp for schedule(static)
      loop_sz_3(sz, i, j, k)
        lr = rop(lr, expr.at(i, j, k));
    } else if (Dim == 4) {
#pragma omp for schedule(static)
      loop_sz_4(sz, i, j, k, l)
        lr = rop(lr, expr.at(i, j, k, l));
    }
#pragma omp critical
    r = rop(r, lr);
  } // end of #omp parallel
  return r;
}

template <class TExpr>
inline auto sum(const TensorExpr<TExpr> &expr) {
  return reduce(expr, std::plus<>(), 0.0);
}

template <class TExpr>
inline auto max(const TensorExpr<TExpr> &expr) {
  using T = std::decay_t<decltype(expr.at(0))>;
  return reduce(expr, NS4_Ops::max_of<>(), std::numeric_limits<T>::lowest());
}

// Important note : this norm flattens the tensor into a euclidean vector
template <class TExpr>
inline auto norm(const TensorExpr<TExpr> &rhs, int nt = 2)
{
  using T = std::decay_t<decltype(rhs.at(0))>;
  T r;
  if(nt == 2)
    r = sqrt(sum(rhs * rhs));
  if(nt == 1)
    r = sum(abs(rhs));
  if(nt == 0)
    r = max(abs(rhs));
  return r;
}

#endif //TENSOREXPR_H
