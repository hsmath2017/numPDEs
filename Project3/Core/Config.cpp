#include <fstream>
#include "Config.h"

#ifndef NDEBUG
int _dbglevel = 0;
std::ofstream _tmpos("tmpos.dat", std::ios::binary);
std::ostream tmpos(_tmpos.rdbuf());
#endif // NDEBUG
