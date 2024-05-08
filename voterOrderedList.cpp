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
    double t(100000);  
    double deltaT(1000); 
    double alpha(0.3);
    int const N1((int)(t/deltaT));
    int const Niter(1);
  
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
    int *rhoo = new int[N1];
    int *survive = new int[N1];
    vector<double> list;
    vector<int> listIndex;
    vector<double> vide;
    vector<int> videIndex;
    double a;

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
    
    for (int k = 0; k < N1; k++)
    {
        rhoo[k] = 0;
        survive[k] = 0;
    } 

    rhoo[0] = Niter*rho_0*(N-rho_0);  
    survive[0] = Niter;    
    
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
		    if (Time[k] < t)
		    {
		    	list.push_back(Time[k]);
	    		listIndex.push_back(k);
		    }    	
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
		
		
		//while (rho != 0 && rho != N)
		while (list.size() > 0)
		{            
		    //cout << list[0] << endl;
		    for (int i = tau/deltaT + 1; i <= list[list.size()-1]/deltaT; i++)
		    {
		    	rhoo[i] += rho*(N-rho);
		    	//cout << i << endl;
		    }
		    
		    if (rho != 0 && rho != N)
		    {
		    	for (int i = tau/deltaT + 1; i <= list[list.size()-1]/deltaT; i++)
				{
					survive[i] += 1;
				}
		    }
		    
		    tau = list[list.size()-1];		   	         
		    
		    Nguy = mtrand1.randInt(N-1);	
		    		    
		    //rho += opinion[Nguy]- opinion[listIndex[list.size()-1]];              
		    //opinion[listIndex[list.size()-1]] = opinion[Nguy];
		  
		    rho += opinion[listIndex[list.size()-1]]-opinion[Nguy];
		    opinion[Nguy] = opinion[listIndex[list.size()-1]];
		
 		    r = mtrand1.rand();
		    list[list.size()-1] += pow(1-r,-1/alpha)-1;
		    //list[list.size()-1] += -log(1-r); 
		    
		    u = list.size()-1;
		   
		    if (list[list.size()-1] > t)
		    {
		    	list.pop_back();
		    	u = 0;
		    }	    
		    
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
    }
    
    ofstream files ("VoterOneRun.dat");
    for (int k = 0; k < N1; k++)
    {
    	files << k*deltaT << " " << 2*rhoo[k]/((double)N*Niter*(N-1)) << " " << survive[k]/(double)Niter << endl;
    }
    
    files.close();
    
    delete[] opinion;
    delete[] Time; 
    delete[] survive;
    delete[] rhoo;  
    delete[] opinion_0;   
}


