#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include "MersenneTwister.h"
using namespace std;



int main()
{  
    MTRand mtrand1;   

    int const N(100);   //number of vertices 
    double alpha(0.3);
    int const Niter(100000);
  
    int j;
    int i(1);
    int u;
    int Nguy;
    int rho_0(0);
    double r;
    int rho(0);
    double *Time = new double[N];
    int *opinion = new int[N];
    int *opinion_0 = new int[N];
    vector<int> rhoo_n(N+1,0);
    vector<double> list;
    vector<int> listIndex;
    vector<double> vide;
    vector<int> videIndex;
    double a;
    vector<int> seen(N,0);
    vector<int> seen_0(N,0);
    vector<int> seen2(N,0);
    vector<int> seen_avg(N,0);

    for (int k = 0; k < N; k++)
    {        
        opinion_0[k] = 0;
        r = mtrand1.rand();
        if (r < 0.5)
        {
            opinion_0[k] = 1;
        }
    }

    for (int k = 0; k < N; k++)
    {
        rho_0 += opinion_0[k];
    }
    
    for (int i = 1; i <= Niter; i++)
    {
    	listIndex = videIndex;
    	list = vide;
    	rho = rho_0;
    	
    	for (int k = 0; k < N; k++)
		{        
		    r = mtrand1.rand();
		    Time[k] = pow(1-r,-1/alpha)-1;
		    //Time[k] = -log(1-r); 	   
		    list.push_back(Time[k]);
	    	listIndex.push_back(k);		        	
	    	opinion[k] = opinion_0[k];   		    		              
		}
		
		for (int k = 2; k <= list.size(); k++)
		{   
		    u = list.size() - k;
		    while (u < list.size()-1 && list[u] < list[u+1])
		    {            
		        a = list[u];  
		        j = listIndex[u];              
		        list[u] = list[u+1]; 
		        listIndex[u] = listIndex[u+1];               
		        list[u+1] = a;
		        listIndex[u+1] = j;
		        u += 1;                           
		    }            
		}
		
		double tau = 0;
		int nbActivation = 0;
		int update = 0;
		int portion = 0;
		seen = seen_0;
		
		while (nbActivation < N)
		{          		    
		    tau = list[list.size()-1];
		    nbActivation++;   		         
		    
		    Nguy = mtrand1.randInt(N-1);
		    
			
										
					    
		    rho += opinion[listIndex[list.size()-1]]-opinion[Nguy];   ////////// moran
		    seen[Nguy]++;
		    seen_avg[Nguy]++;
		    opinion[Nguy] = opinion[listIndex[list.size()-1]];		

/*		    		    
		    rho += opinion[Nguy]- opinion[listIndex[list.size()-1]];   ////////////  voter 
		    seen[listIndex[list.size()-1]]++; 
		    seen_avg[listIndex[list.size()-1]] += 1;          
		    opinion[listIndex[list.size()-1]] = opinion[Nguy];
*/		    
		    
 		    r = mtrand1.rand();
		    list[list.size()-1] += pow(1-r,-1/alpha)-1;
		    //list[list.size()-1] += -log(1-r); 
		    
		    u = list.size()-1;   
		    
		    while (list[u] > list[u-1] && u > 0)
		    {
		    	a = list[u];  
		        j = listIndex[u];              
		        list[u] = list[u-1]; 
		        listIndex[u] = listIndex[u-1];               
		        list[u-1] = a;
		        listIndex[u-1] = j;
		        u -= 1;
		    } 
		}
		
		for (int k = 0; k < N; k++)
		{
			seen2[k] += pow(seen[k]-1,2); 
		}		
	}
	
		
    
    ofstream fil ("voterNumberOfShots.dat");
    
	for (int i = 0; i < N; i++)
	{
	    fil << i << " " << seen_avg[i]/(double)Niter << " " << seen2[i]/((double)Niter) << endl;
	}
    
    fil.close();
    
    delete[] opinion;
    delete[] Time; 
    delete[] opinion_0;   
}


