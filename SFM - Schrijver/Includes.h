
#ifndef Includes_h
#define Includes_h

#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <numeric>
#include <unordered_set>

using namespace std;

typedef unordered_set<int> Subset;
typedef vector<int> VI;
typedef vector<double> vec;

#define PB push_back
#define MP make_pair
#define F first
#define S second
#define SZ(a) (int)(a.size())
#define SET(a,b) memset(a,b,sizeof(a))
#define LET(x,a) __typeof(a) x(a)
#define ALL(c) (c).begin(),(c).end()

#define TRACE
#ifdef TRACE
#define trace(...) __f(#__VA_ARGS__, __VA_ARGS__)
template <typename Arg1>
void __f(const char* name, Arg1&& arg1)	{
    cerr << name << " : " << arg1 << std::endl;}
template <typename Arg1, typename... Args>
void __f(const char* names, Arg1&& arg1, Args&&... args){
    const char* comma = strchr(names + 1, ',');cerr.write(names, comma - names) << " : " << arg1<<" | ";__f(comma+1, args...);}
#else
#define trace(...)
#endif



#endif /* Includes_h */
