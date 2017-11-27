
#ifndef Schrijver_hpp
#define Schrijver_hpp

#include "Includes.h"
#include "Graph.hpp"

class SF
{
public:
    
    /*
     - change definition of evalute function as per requirement/problem
     */
    
    VI modular;
    int numNodes;
    VI V;   // ground set
    SF(VI mod_init); // change arg. list as per problem requirement
    double evaluate(Subset s);
    Subset minimize();
};


class totalOrder
{
public:
    VI order;   // total order
    vec greedyVector;   // corresponding greedy vector
    
    totalOrder(VI init,SF problem);
    void display();
    int gapST(int s,int t);
    vector<totalOrder> generateNewOrder(int s,int t,SF problem);
    vec arrangeGreedyVector(vec v); // arranges v , ordered as increasing order of vertices,  as per the order
};


#endif /* Schrijver_hpp */
