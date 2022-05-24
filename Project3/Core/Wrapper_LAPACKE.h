#ifndef WRAPPER_LAPACKE_H
#define WRAPPER_LAPACKE_H

#ifdef USE_MKL
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#else
#ifdef USE_AOCL
#include <cblas.h>
#include <lapacke.h>
#else
#ifdef USE_OPENBLAS
#include <cblas.h>
#include <lapacke.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif // USE_OPENBLAS
#endif // USE_AOCL
#endif // USE_MKL


#endif //WRAPPER_LAPACKE_H
