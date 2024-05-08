#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "MersenneTwister.h"
#include <vector>
//#include <boost/math/special_functions/erf.hpp>


using namespace std;

int main()
{   
    MTRand mtrand1;

    long long int const N(100);
    int Niter = 1000000;
    double alpha = 0.7;
    
    double r;
    double t = 100000;
    double deltaT = 1;
    int probamax = t/deltaT;    
    vector<double> activationTime(N,0);
	vector<double> lastInteractionTime(N,0);
	vector<double> act(N);
	vector<double> list;
	vector<int> qui;

	double timeMin;
	int collegue;
	int indicateur = 1;
	
    
    for (int u = 0; u < N; u++)
	{
		act[u] = 1;
		r = mtrand1.randDblExc();
		//activationTime[u] = (pow(1-r,-1/alpha)-1)/act[u];
	}
	
	 
	for (int u = 0; u < Niter; u++)
	{
		r = mtrand1.randDblExc();
		activationTime[0] += (pow(1-r,-1/alpha)-1)/act[0];
		list.push_back(activationTime[0]);
	}	
	
	for (int i = 1; i < N; i++)
	{			
		while (activationTime[i] < activationTime[0])
		{
			r = mtrand1.randDblExc();
			activationTime[i] += (pow(1-r,-1/alpha)-1)/act[i];
		
			if (mtrand1.randInt(N-2) == 0)
			{
				list.push_back(activationTime[i]);
			}
  			//list.push_back(activationTime[i]);
		}
	}
/*	
	double maximus = 0;
	
	for (int k = 0; k < list.size() ; k++)
    {
    	if (list[k] > maximus)
    	{
    		maximus = list[k];
    	}
    }
*/    
  
	vector<double> proba(1000000,0);
		
	
	cout << "début du tri de la liste qui fait " << list.size()  << endl;
	
	// ordonner la liste
	int u,a,j;
	for (int k = 1; k < list.size() ; k++)
    {   
        u = k;
        while (u > 0 && list[u] < list[u-1])
        {            
            a = list[u];  
            //j = qui[u];              
            list[u] = list[u-1]; 
            //qui[u] = qui[u-1];               
            list[u-1] = a;
            //qui[u-1] = j;
            u -= 1;                           
        } 
        //cout << k << endl;           
    }
    
	int normal = 0;
	int integer;
	for (int j = 1; j < list.size(); j++)
	{
		integer = (int)((list[j]-list[j-1])/deltaT);
		if (integer < proba.size())
		{
			proba[integer] += 1;
			normal++;	
		}
		
	}		
    	
	ofstream filo ("intereventDistr.dat");
    for (int i = 1; i < proba.size(); i++)
    {
    	filo << deltaT*i << " " << proba[i]/(normal*deltaT) << endl;
    }
    filo.close();   
}



