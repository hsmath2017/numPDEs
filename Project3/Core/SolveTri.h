#ifndef _SOLVETRI_H_
#define _SOLVETRI_H_

#include "Core/Config.h"
#include <vector>
#include <assert.h>

std::vector<Real> solveL(const std::vector<Real> &a, const std::vector<Real> &p, const std::vector<Real> &rhs)
{
    assert((int)a.size() + 1 == (int)p.size() && (int)p.size() == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[0] = rhs[0] / p[0];
    for (int i = 1; i < N; i++)
    {
        result[i] = (rhs[i] - result[i - 1] * a[i - 1]) / p[i];
    }
    return result;
}

std::vector<Real> solveU(const std::vector<Real> &q, const std::vector<Real> &rhs)
{
    assert((int)q.size() + 2 == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[N - 1] = rhs[N - 1];
    result[N - 2] = rhs[N - 2];
    for (int i = N - 3; i >= 0; i--)
    {
        result[i] = rhs[i] - q[i] * result[i + 1];
    }
    return result;
}

std::vector<Real> solveD(const std::vector<Real> &t, Real s, const std::vector<Real> &rhs)
{
    assert((int)t.size() + 1 == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[N - 1] = (rhs[N - 1] - rhs[0] * s) / (1 - s * t[0]);
    for (int i = N - 2; i >= 0; i--)
    {
        result[i] = rhs[i] - t[i] * result[N - 1];
    }
    return result;
}

std::vector<Real> solvePeriodicTri(const std::vector<Real> &a, const std::vector<Real> &b, const std::vector<Real> &c, const std::vector<Real> &rhs)
{
    assert(a.size() == b.size() && b.size() == c.size() && c.size() == rhs.size());
    int N = rhs.size();

    //get a in L
    std::vector<Real> newa;
    newa.assign(a.begin() + 1, a.end());

    //get p, q , r , s in L and Utilde
    std::vector<Real> p(N);
    std::vector<Real> q(N - 2);
    std::vector<Real> r(N - 1);
    p[0] = b[0];
    q[0] = c[0] / p[0];
    r[0] = a[0] / p[0];
    for (int i = 1; i < N - 2; i++)
    {
        p[i] = b[i] - a[i] * q[i - 1];
        q[i] = c[i] / p[i];
        r[i] = -a[i] * r[i - 1] / p[i];
    }
    p[N - 2] = b[N - 2] - a[N - 2] * q[N - 3];
    r[N - 2] = (c[N - 2] - a[N - 2] * r[N - 3]) / p[N - 2];
    p[N - 1] = b[N - 1] - a[N - 1] * r[N - 2];
    Real s = c[N - 1] / p[N - 1];

    //get t in D
    std::vector<Real> t(N - 1);
    t[N - 2] = r[N - 2];
    for (int i = N - 3; i >= 0; i--)
    {
        t[i] = r[i] - q[i] * t[i + 1];
    }

    std::vector<Real> tmp;
    tmp = solveL(newa, p, rhs);
    tmp = solveU(q, tmp);
    return solveD(t, s, tmp);
}

std::vector<Real> solveLowerTri(const std::vector<Real> &a, const
                           std::vector<Real> &p, const std::vector<Real> &rhs)
{
    assert((int)a.size() + 1 == (int)p.size() && (int)p.size() == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[0] = rhs[0] / p[0];
    for (int i = 1; i < N; i++)
    {
        result[i] = (rhs[i] - result[i - 1] * a[i - 1]) / p[i];
    }
    return result;
}
 
std::vector<Real> solveUpperTri(const std::vector<Real> &q, const std::vector<Real> &rhs)
{
    assert((int)q.size() + 1 == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[N - 1] = rhs[N - 1];
    for (int i = N - 2; i >= 0; i--)
    {
        result[i] = rhs[i] - q[i] * result[i + 1];
    }
    return result;
}
 
 
std::vector<Real> solveLsp(const std::vector<Real> &a, const
                           std::vector<Real> &p, const Real e, const std::vector<Real> &rhs)
{
    assert((int)a.size() + 1 == (int)p.size() && (int)p.size() == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[0] = rhs[0] / p[0];
    for (int i = 1; i < N - 1; i++)
    {
        result[i] = (rhs[i] - result[i - 1] * a[i - 1]) / p[i];
    }
    result[N - 1] = (rhs[N - 1] - result[N - 2]*a[N - 2] - result[N -
                           3]*e) / p[N - 1];
    return result;
}
 
std::vector<Real> solveUsp(const std::vector<Real> &q, const Real d, const std::vector<Real> &rhs)
{
    assert((int)q.size() + 1 == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[N - 1] = rhs[N - 1];
    for (int i = N - 2; i >= 1; i--)
    {
        result[i] = rhs[i] - q[i] * result[i + 1];
    }
    result[0] = rhs[0] - q[0] * result[1] - d * result[2];
    return result;
}
 
std::vector<Real> solveTri(const std::vector<Real> &a, const
                           std::vector<Real> &b, const
                           std::vector<Real> &c, const std::vector<Real> &rhs)
{
     assert(b.size() == rhs.size() && a.size() == b.size() - 1 &&
            c.size() == b.size() - 1);
     int N = rhs.size();
     
     //get a in L
     std::vector<Real> newa;
     newa.assign(a.begin(), a.end());
 
     //get p, q in L and U
     std::vector<Real> p(N);
     std::vector<Real> q(N - 1);
     p[0] = b[0];
     q[0] = c[0] / p[0];
     for (int i = 1; i < N - 1; i++)
     {
          p[i] = b[i] - a[i - 1] * q[i - 1];
          q[i] = c[i] / p[i];
     }
     p[N - 1] = b[N - 1] - a[N - 2] * q[N - 2];
 
    std::vector<Real> tmp;
    tmp = solveLowerTri(newa, p, rhs);
    tmp = solveUpperTri(q, tmp);
    return tmp;
}
 
std::vector<Real> solveTrisp(const std::vector<Real> &a, const
                             std::vector<Real> &b, const
                             std::vector<Real> &c, const Real d, const
                             Real e, const std::vector<Real> &rhs)
{
     assert(b.size() == rhs.size() && a.size() == b.size() - 1 &&
            c.size() == b.size() - 1);
     int N = rhs.size();
     
     //get a in L, with a[N-2] to be modify
     std::vector<Real> newa;
     newa.assign(a.begin(), a.end());
 
     //get p, q ,d ,e in L and U
     std::vector<Real> p(N);
     std::vector<Real> q(N - 1);
     p[0] = b[0];
     q[0] = c[0] / p[0];
     Real newd = d / p[0];
     p[1] = b[1] - a[0] * q[0];
     q[1] = (c[1] - newd * a[0]) / p[1];
     for (int i = 2; i < N - 1; i++)
     {
          p[i] = b[i] - a[i - 1] * q[i - 1];
          q[i] = c[i] / p[i];
     }
     Real newe = e;
     newa[N - 2] = a[N - 2] - newe * q[N - 3];
     p[N - 1] = b[N - 1] - newa[N - 2] * q[N - 2];
 
    std::vector<Real> tmp;
    tmp = solveLsp(newa, p, newe, rhs);
    tmp = solveUsp(q, newd, tmp);
    return tmp;
}

#endif
