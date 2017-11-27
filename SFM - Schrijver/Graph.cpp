
#include "Graph.hpp"


Graph::Graph(int n)
{
    numNodes = n;
    g.resize(n);
}

Graph::Graph(int n,vector<VI> orderings)
{
    numNodes = n;
    g.resize(n);
    
    for(int i=0;i<SZ(orderings);i++)
    {
        for(int j=0;j<SZ(orderings[i]);j++)
        {
            for(int k=j+1;k<SZ(orderings[i]);k++)
            {
                g[orderings[i][j]].PB(orderings[i][k]);
            }
        }
    }
    
}

bool Graph::existPath(Subset P,Subset N)
{
    queue<int> q;
    vector<bool> vis(numNodes+1,0);
    
    for(auto w:P)
        q.push(w),vis[w]=1;
    
    // using multi-source bfs to find existance of path from P to N
    
    while(!q.empty())
    {
        int u = q.front();
        q.pop();
        
        for(auto w:g[u])
        {
            if(!vis[w])
            {
                if(N.find(w)!=N.end())
                    return true;
                vis[w]=1;
                q.push(w);
            }
        }
    }
    return false;
}

void Graph::distance(Subset P, vector<int> &d)
{
    d.resize(numNodes+1);
    fill(ALL(d),INT_MAX);
    queue<int> q;
    vector<bool> vis(numNodes+1,0);
    
    for(auto w:P)
        q.push(w),d[w]=0,vis[w]=1;
    
    while(!q.empty())
    {
        int u = q.front();
        q.pop();
        
        for(auto w:g[u])
        {
            if(!vis[w])
            {
                d[w]=min(d[w],d[u]+1);
                vis[w]=1;
                q.push(w);
            }
        }
    }
    
}

bool Graph::existEdge(int s,int t)
{
    for(auto w:g[s])
    {
        if(w==t)
            return true;
    }
    return false;
}
