
#include "Schrijver.hpp"
#include "Graph.hpp"
#include "Reduce.hpp"

SF::SF(VI mod_init)
{
    modular = mod_init;
    numNodes = SZ(modular);
    for(int i=0;i<SZ(modular);i++)
        V.PB(i);
}

double SF::evaluate(Subset s)
{
    double ans = (SZ(V)-SZ(s))*SZ(s);
    
    for(auto w:s)
        ans+=modular[w];
    
    return ans;
}

totalOrder::totalOrder(VI init,SF problem): order(init)
{
    Subset s;
    greedyVector.resize(SZ(init));
    double current,previous = problem.evaluate(s);

    for(auto w:order)
    {
        s.insert(w);
        current = problem.evaluate(s);
        greedyVector[w] = current - previous;
        previous = current;
    }
    
}
void totalOrder::display()
{
    cout<<"The total ordering is:"<<endl;
    cout<<"[";
    for(auto w:order)
        cout<<w<<", ";
    cout<<"]"<<endl;
    
    cout<<"The corresponding greedy vector is:"<<endl;
    cout<<"[";
    for(auto w:greedyVector)
        cout<<w<<", ";
    cout<<"]"<<endl;
    
}

int totalOrder::gapST(int s,int t)
{
    int s_pos = (int)(find(ALL(order),s) - order.begin());
    int t_pos = (int)(find(ALL(order),t) - order.begin());
    
    if(s_pos <= t_pos)
    {
        return t_pos-s_pos;
    }
    return 0;
    
}

/*
 - helper function to generate all total orders
 - after shifting u , contained in (s,t] , to 
 - just before s
 */
vector<totalOrder> totalOrder::generateNewOrder(int s, int t,SF problem)
{
    int s_pos = (int)(find(ALL(order),s) - order.begin());
    int t_pos = (int)(find(ALL(order),t) - order.begin());
    
    vector<totalOrder> answerOrders;
    
    for(int i=s_pos+1;i<=t_pos;i++)
    {
        VI tmpOrder(SZ(order),0);
        
        int it1=0,it2=0;
        while(it1<SZ(order) && it2<s_pos)
        {
            tmpOrder[it1++]=order[it2++];
        }
        tmpOrder[it1++]=order[i];
        while(it1<SZ(order) && it2<SZ(order))
        {
            if(it2==i)
                it2++;
            else
            {
                tmpOrder[it1++]=order[it2++];
            }
        }
        answerOrders.PB(totalOrder(tmpOrder,problem));
    }
    return answerOrders;
}

/*
 - helper function to arrange the greedy vector
 - in the order specified by the total order
 -
 - by default, the greedy vector is in increasing order of elements i.e. [0,1,2,3,....|V|]
 - to maintain uniformity
 */
vec totalOrder::arrangeGreedyVector(vec v)
{
    vec answer(SZ(v));
    
    for(int i=0;i<SZ(order);i++)
    {
        answer[i]=v[order[i]];
    }
    return answer;
}


// helper function to create first dummy ordering
VI makeDummyOrder(int n)
{
    VI v;
    for(int i=0;i<n;i++)
        v.PB(i);
    random_shuffle(ALL(v));
    
    return v;
}

/*
 - main function to pick the subset corresponding to minimization of SF
 */
Subset SF::minimize()
{
    Subset picked;
    
    // make the needed variables
    vec x,new_x;
    Subset P;
    Subset N;
    Subset Neutral;
    vector<int> d;
    vector<totalOrder> B;
    vector<double> W;
    Graph  orderingGraph;
    int dt=INT_MAX;
    int t=INT_MAX;
    int s=INT_MAX;
    int alpha=INT_MAX;
    int beta=INT_MAX;
    
    // make a initial dummy ordering
    B.PB(totalOrder(makeDummyOrder(SZ(V)),*this));
    x = B[0].greedyVector;
    W.PB(1);
    
    while(true)
    {
        
        //constructing graph
        orderingGraph  = Graph(numNodes);
        for(auto w:B)
        {
            for(int i=0;i<SZ(w.order);i++)
            {
                for(int j=i+1;j<SZ(w.order);j++)
                {
                    orderingGraph.g[w.order[i]].PB(w.order[j]);
                }
            }
        }
        
        // finding P,N, Neutral (enteries with 0)
        P.clear();N.clear();Neutral.clear();
        for(int i=0;i<SZ(x);i++)
        {
            if(x[i]<0)
                N.insert(i);
            else if(x[i]>0)
                P.insert(i);
            else
                Neutral.insert(i);
        }
        
        // not path from P to N, breaking condition
        // find all vertices which have a path to N
        if(!orderingGraph.existPath(P,N))
        {
            picked = N;
            for(auto w:Neutral)
            {
                Subset tmp; tmp.insert(w);
                if(orderingGraph.existPath(tmp, N))
                    picked.insert(w);
            }
            return picked;
        }
        
        d.clear();
        orderingGraph.distance(P, d);
        dt=0;
        
        //finding t
        for(int i=0;i<numNodes;i++)
        {
            if(N.find(i)!=N.end() && d[i]!=INT_MAX && d[i]>=dt)
            {
                dt=d[i];
                t=i;
            }
        }
        
        //finding s
        for(int i=0;i<numNodes;i++)
        {
            if(orderingGraph.existEdge(i, t) && d[i]==(d[t]-1) )
                s=i;
        }
        
        // finding the total order with maximum separation between s and t
        int mxx=0;
        int idxB = 0;   // index of the selected total order in this iteration
        alpha=0;beta=0;
        
        for(int i=0;i<SZ(B);i++)
        {
            int tmp = B[i].gapST(s,t);
            if(tmp > mxx)
            {
                mxx = tmp;
                alpha = mxx;
                idxB = i;
                beta = 1;
            }
            else if(tmp == mxx)
            {
                beta++;
            }
        }
        
        // make all new totalOrders from the selected totalOrdering
        double delta = 0.0;
        vector<totalOrder> new_B = B[idxB].generateNewOrder(s, t, *this);
        vector<double> new_W(SZ(new_B),0);
        totalOrder &base_B = B[idxB];
        vec ref = base_B.arrangeGreedyVector(base_B.greedyVector);
        int s_pos = int(find(ALL(base_B.order),s)-base_B.order.begin());
        int t_pos = int(find(ALL(base_B.order),t)-base_B.order.begin());
        
        vector<vec> M(SZ(new_B));
        bool flagZeroDelta = false;
        int  idxZeroDelta = 0;
        
        for(int  i=0;i<SZ(new_B);i++)
        {
            
            M[i] = base_B.arrangeGreedyVector(new_B[i].greedyVector);
            for(int j=0;j<SZ(M[i]);j++)
            {
                M[i][j] = M[i][j] - ref[j];
                if(j==(s_pos+1+i) && M[i][j]<numeric_limits<double>::epsilon())
                {
                    flagZeroDelta = true;
                    idxZeroDelta = i;
                }
            }
        }
        
        if(flagZeroDelta)
        {
            delta = 0;
            new_W[idxZeroDelta]=1;
        }
        else
        {
            for(int i=(SZ(M)-1);i>=0;i--)
            {
                double rhs=0.0;
                if(i==(SZ(M)-1))
                    rhs=1.0;
                
                int j = s_pos + (i+1);
                double lhs=0.0;
                for(int k=i+1;k<SZ(M);k++)
                {
                    lhs += M[k][j]*new_W[k];
                }
                new_W[i] = (rhs - lhs)/M[i][j];
            }
            delta = 0.0;
            for(auto w:new_W)
                delta +=w;
            delta = 1.0/delta;
            for(auto &w:new_W)
                w *= delta*w;
        }// end of else
            // find y and x'
            vec y=x;
            double lambda = W[idxB];
            double multiplier = 1.0;
            y[s] -= lambda*delta;
            y[t] += lambda*delta;
            
            if (y[t] <= 0)
            {
                new_x = y;
                W[idxB] = 0;
            }
            else
            {
                new_x = x;
                new_x[s] += new_x[t];
                multiplier = -new_x[t] / (delta*lambda);
                new_x[t] = 0;
                W[idxB] = lambda * (1-multiplier);
            }
            
            
            // Add the new ordering and the related weights to the maintained list
            for (totalOrder w: new_B)
                B.PB(w);
    
            for (double w: new_W) {
                W.PB(w*lambda*multiplier);
            }
        
        
            //remove total orders with zero weight
            int pos=0;
            while(pos<SZ(B))
            {
                if(W[pos]==0.0)
                {
                    B.erase(B.begin() + pos);
                    W.erase(W.begin() + pos);
                }
                else
                    pos++;
            }
            
            // remove all affine dependent total orders
            
            bool ok = false;
            do
            {
                ok = false;
                vector<vec> tmp;
                for(auto w:B)
                    tmp.PB(w.greedyVector);
                
                reduce(tmp,W);
                
                int pos=0;
                while(pos<SZ(B))
                {
                    if(W[pos]==0.0)
                    {
                        B.erase(B.begin() + pos);
                        W.erase(W.begin() + pos);
                        ok=true;
                    }
                    else
                        pos++;
                }
                
            }while(ok);
            
            x= new_x;
        
        
    }// end of iteration loop !!
    
    return picked;
}

