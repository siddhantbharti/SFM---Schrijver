
#include "Schrijver.hpp"



/*
 - it is a brute O(2^n) function 
 - to check correctness
 */
pair<Subset,double> solve(int n, VI modular)
{
    double ansS=numeric_limits<double>::max();
    Subset ansF;
    
    for(int mask=0;mask<(1LL<<n);mask++)
    {
        int n1=0;
        double tmpS=0.0;
        Subset tmpF;
        for(int i=0;i<n;i++)
        {
            if(mask&(1LL<<i))
                n1++,tmpS+=modular[i],tmpF.insert(i);
        }
        tmpS+=(n-n1)*n1;
        if(tmpS<ansS)
        {
            ansS=tmpS;
            ansF=tmpF;
        }
    }
    
    return {ansF,ansS};
}

/*
 - generates modular vector to be added (corresponding index elements) to submodular function
 */
VI generateModular(int n)
{
    VI ans;
    for(int i=1;i<=n;i++)
        ans.PB((random()%100)*(((i&1)==0)?1:-1));
    random_shuffle(ALL(ans));
    return ans;
}


int main(int argc, const char * argv[])
{
    // n is the size of ground set to be considered
    
    int n;
    cin>>n;
    Subset ans;
    
    VI modular = generateModular(n);
    
    cout<<"Schrijver: "<<endl;
    SF problem(modular);
    ans = problem.minimize();
    for(auto w:ans)
        cout<<w<<" ";
    cout<<endl;
    cout<<problem.evaluate(ans)<<endl;
    
    cout<<"\nBrute Force: "<<endl;
    pair<Subset,int> ans1 = solve(n,modular);
    for(auto w:ans1.F)
        cout<<w<<" ";
    cout<<endl;
    cout<<ans1.S<<endl;
	return 0;
}




