#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>

using Real = double;

#ifndef DIM
#define DIM 2
#endif

#if DIM == 2
const int SpaceDim = 2;
#elif DIM == 3
const int SpaceDim = 3;
#endif

//=================================================

// debug output issues
#ifndef NDEBUG

extern int _dbglevel;
extern std::ostream tmpos;
#define push_dbglevel(x)  _dbglevel = (_dbglevel << 4) | (x)
#define pop_dbglevel()    _dbglevel = _dbglevel >> 4
#define reset_dbglevel(x) _dbglevel = (x)
#define get_dbglevel()    (_dbglevel & 0x0f)

#else

#define tmpos 0 && std::cout
#define push_dbglevel(x)
#define pop_dbglevel()
#define reset_dbglevel(x)
#ifndef DBGLEVEL
#define DBGLEVEL -1
#endif
#define get_dbglevel() (DBGLEVEL)

#endif // NDEBUG

#define dbgcout  (get_dbglevel()>=0) && std::cout
#define dbgcout1 (get_dbglevel()>=1) && std::cout << "  "
#define dbgcout2 (get_dbglevel()>=2) && std::cout << "    "
#define dbgcout3 (get_dbglevel()>=3) && std::cout << "      "

#endif // CONFIG_H
