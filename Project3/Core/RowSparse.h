#ifndef ROWSPARSE_H
#define ROWSPARSE_H

#include <vector>
#include <initializer_list>
#include <functional>
#include <algorithm>
#include "Core/Tensor.h"

template <class T_RowIndex = int, class T_ColIndex = T_RowIndex>
class RowSparse
{
public:
  RowSparse() = default;

  std::size_t insertRow(const T_RowIndex &r, std::initializer_list<T_ColIndex> cols, std::initializer_list<Real> vals);

  std::size_t insertRow(const T_RowIndex &r, int nz, const T_ColIndex *cols, const Real *vals);

  std::size_t getNumRows() const { return rows.size(); }

  // calculate z = alpha * Ax + y
  template <class T1, class T2, class T3>
  void AXPY(Real alpha, const T1 &x, const T2 &y, T3 &z) const;

  // calculate z = alpha * ax + y, compressed row output
  template <class T1, class T2>
  Tensor<Real,1> AXPY(Real alpha, const T1 &x, const T2 &y) const;

  // iterators
public:
  class RowSparseIterator {
    typename std::vector<Real>::const_iterator       pVal;
    typename std::vector<T_ColIndex>::const_iterator pCol;
    typename std::vector<int>::const_iterator        pNz;
    typename std::vector<T_RowIndex>::const_iterator pRw;

  public:
    friend class     RowSparse;
    const Real       &getValue(int j)  const { return *(pVal+j); }
    const T_ColIndex &getColumn(int j) const { return *(pCol+j); }
    const T_RowIndex &getRow()         const { return *pRw; }
    int              getNz()           const { return *pNz; }

    RowSparseIterator &operator++() {
      pVal += getNz();
      pCol += getNz();
      ++pNz;
      ++pRw;
      return *this;
    }
    void operator++(int) { ++(*this); }
    bool operator!=(const RowSparseIterator &rhs) const { return this->pRw != rhs.pRw; }
    bool operator==(const RowSparseIterator &rhs) const { return this->pRw == rhs.pRw; }
    std::ptrdiff_t operator-(const RowSparseIterator &rhs) const { return this->pRw - rhs.pRw; }
  };

  RowSparseIterator cbegin() const {
    RowSparseIterator it;
    it.pVal = values.cbegin();
    it.pCol = columns.cbegin();
    it.pNz = nonZeros.cbegin();
    it.pRw = rows.cbegin();
    return it;
  }

  RowSparseIterator cend() const {
    RowSparseIterator it;
    it.pRw = rows.cend();
    return it;
  }

  template <class T_Less = std::less<T_RowIndex>>
  RowSparseIterator find_fast(const T_RowIndex &ri, const T_Less &ls) const;

protected:
  std::vector<Real>       values;
  std::vector<T_ColIndex> columns;
  std::vector<int>        nonZeros;
  std::vector<int>        rowBegin;
  std::vector<T_RowIndex> rows;
};

//============================================================

template <class T_RowIndex, class T_ColIndex>
inline
std::size_t RowSparse<T_RowIndex, T_ColIndex>::insertRow(const T_RowIndex &r, int nz, const T_ColIndex *cols, const Real *vals)
{
  std::size_t a = rows.size();
  rows.push_back(r);
  rowBegin.push_back(values.size());
  nonZeros.push_back(nz);
  values.insert(values.cend(), vals, vals + nz);
  columns.insert(columns.cend(), cols, cols + nz);
  return a;
}

template <class T_RowIndex, class T_ColIndex>
inline
std::size_t RowSparse<T_RowIndex, T_ColIndex>::insertRow(const T_RowIndex &r,
                                                         std::initializer_list<T_ColIndex> cols,
                                                         std::initializer_list<Real> vals)
{
  return insertRow(r, cols.size(), cols.begin(), vals.begin());
}

template <class T_RowIndex, class T_ColIndex>
template <class T1, class T2, class T3>
inline
void RowSparse<T_RowIndex, T_ColIndex>::AXPY(Real alpha, const T1 &x, const T2 &y, T3 &z) const
{
  for(auto rowIt = cbegin(); rowIt != cend(); ++rowIt) {
    Real Ax = 0.0;
    for(int k=0; k<rowIt.getNz(); ++k)
      Ax += rowIt.getValue(k) * x(rowIt.getColumn(k));
    z(rowIt.getRow()) = alpha * Ax + y(rowIt.getRow());
  }
}

template <class T_RowIndex, class T_ColIndex>
template <class T1, class T2>
inline
Tensor<Real,1> RowSparse<T_RowIndex, T_ColIndex>::AXPY(Real alpha, const T1 &x, const T2 &y) const
{
  Tensor<Real,1> axpy(getNumRows());
  int i = 0;
  for(auto rowIt = cbegin(); rowIt != cend(); ++rowIt, ++i) {
    Real Ax = 0.0;
    for(int k=0; k<rowIt.getNz(); ++k)
      Ax += rowIt.getValue(k) * x(rowIt.getColumn(k));
    axpy(i) = alpha * Ax + y(rowIt.getRow());
  }
  return axpy;
}

template <class T_RowIndex, class T_ColIndex>
template <class T_Less>
inline
auto RowSparse<T_RowIndex, T_ColIndex>::find_fast(const T_RowIndex &ri, const T_Less &ls) const -> RowSparseIterator
{
  auto it = std::lower_bound(rows.cbegin(), rows.cend(), ri, ls);
  if(it == rows.cend() || ls(ri, *it))
    return cend();
  auto i = it - rows.cbegin();
  RowSparseIterator rsit;
  rsit.pRw = it;
  rsit.pNz = nonZeros.cbegin() + i;
  rsit.pCol = columns.cbegin() + rowBegin[i];
  rsit.pVal = values.cbegin() + rowBegin[i];
  return rsit;
}

#endif //ROWSPARSE_H
