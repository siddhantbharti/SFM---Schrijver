
#include "Reduce.hpp"
#include "Includes.h"

T GaussJordan(VVT &a, VVT &b)
{
    const int n=a.size(),m=b[0].size();
    
    VI irow(n),icol(n),ipiv(n);
    T det=1;
    for(int i=0;i<n;i++)
    {
        int pj=-1,pk=-1;
        for(int j=0;j<n;j++)if(!ipiv[j])
            for(int k=0;k<n;k++)if(!ipiv[k])
                if(pj==-1||fabs(a[j][k])>fabs(a[pj][pk]))
                {pj=j;pk=k;}
        
        if(fabs(a[pj][pk])<EPS){return 1;}
        
        ipiv[pk]++;
        swap(a[pj], a[pk]);
        swap(b[pj],b[pk]);
        if(pj!=pk)det*=-1;
        irow[i]=pj;
        icol[i]=pk;
        T c =1.0/a[pk][pk];
        det*=a[pk][pk];a[pk][pk]=1.0;
        
        for(int p=0;p<n;p++)a[pk][p]*=c;
        for(int p=0;p<m;p++)b[pk][p]*=c;
        
        for(int p=0;p<n;p++)if(p!=pk)
        {
            c=a[p][pk];a[p][pk]=0;
            for(int q=0;q<n;q++)a[p][q]-=a[pk][q]*c;
            for(int q=0;q<m;q++)b[p][q]-=b[pk][q]*c;
        }
    }
    for(int p=n-1;p>=0;p--)
        if(irow[p]!=icol[p])
            for(int k=0;k<n;k++)
                swap(a[k][irow[p]],a[k][icol[p]]);
    return det;
}

void reduce(vector<vec> &B,vec &W)
{
    // making working space from Gauss-Jordan Elimination
    VVT a;
    VVT b;
    int n = max(SZ(B[0])+1,SZ(B));
    
    for(int i=0;i<n;i++)
    {
        VT tmp(1,0.0);
        b.PB(tmp);
        tmp.resize(n,0.0);
        a.PB(tmp);
    }
    
    for(int i=0;i<SZ(B);i++)
    {
        for(int j=0;j<SZ(B[0]);j++)
            a[j][i]=B[i][j];
        a[SZ(B[0])][i]=1;
    }
    
    bool exiting = GaussJordan(a, b);
    if(exiting)
        return;
    
    double ratio = numeric_limits<double>::max();
    
    for(int i=0;i<n;i++)
    {
        if(b[i][0]>0)
            ratio = min(ratio,W[i]/b[i][0]);
    }
    if(ratio !=numeric_limits<double>::max())
    {
        for(int i=0;i<n;i++)
        {
            W[i]=W[i]-ratio*b[i][0];
        }
    }
    
}
