
#ifndef Graph_hpp
#define Graph_hpp

#include "Includes.h"

/*
 - using adjacency list implementation of graph
 */

class Graph
{
public:
    int numNodes;
    vector<vector<int> > g;
    
    Graph(){};
    Graph(int n);
    Graph(int n,vector<VI> orderings);
    bool existPath(Subset P,Subset N);
    void distance(Subset P,vector<int> &d);
    bool existEdge(int s,int t);
};


#endif /* Graph_hpp */
