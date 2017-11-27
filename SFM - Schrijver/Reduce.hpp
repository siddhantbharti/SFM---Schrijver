

#ifndef Reduce_hpp
#define Reduce_hpp

#include "Includes.h"

const double EPS = 1e-10;
typedef double T;typedef vector<T> VT;typedef vector<VT> VVT;

void reduce(vector<vec> &B, vec &w);
T GaussJordan(VVT &a, VVT &b);

#endif /* Reduce_hpp */
