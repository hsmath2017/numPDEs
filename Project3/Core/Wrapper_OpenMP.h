#ifndef WRAPPER_OPENMP_H
#define WRAPPER_OPENMP_H

#ifndef OMP_PAR_LOWER_BOUND
#define OMP_PAR_LOWER_BOUND 256
#endif // OMP_PAR_LOWER_BOUND

const int OMP_Par_Lower_Bound = OMP_PAR_LOWER_BOUND;

extern int numInitialThreads;

#ifdef USE_OPENMP

#include <omp.h>

class OMPSwitch
{
public:
  OMPSwitch(int targetThreads, bool e = true) : doSwitch(e)
  {
    if(doSwitch) {
      oldThreads = omp_get_max_threads();
      omp_set_num_threads(targetThreads);
    }
  }

  ~OMPSwitch() {
    if(doSwitch)
      omp_set_num_threads(oldThreads);
  }

protected:
  bool doSwitch;
  int  oldThreads;
};

#else

class OMPSwitch
{
public:
  OMPSwitch(int, bool = true) { }
};

#endif // USE_OPENMP

#endif //WRAPPER_OPENMP_H
