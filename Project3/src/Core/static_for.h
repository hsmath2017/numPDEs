#ifndef STATIC_FOR_H
#define STATIC_FOR_H

#include <type_traits>

template <int initVal, int endVal>
struct static_for {
  template <class Lambda, class Ts>
  static auto execute(Lambda ld, const Ts& arg) {
    return static_for<initVal+1,endVal>::execute(ld, ld(std::integral_constant<int,initVal>(), arg));
  }
};

template <int endVal>
struct static_for<endVal,endVal> {
  template <class Lambda, class Ts>
  static auto execute(Lambda ld, const Ts& arg) { return arg; }
};

#endif //STATIC_FOR_H
