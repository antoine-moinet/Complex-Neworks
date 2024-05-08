#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include "MersenneTwister.h"
using namespace std;



int main()
{  
    MTRand randdd;
	int Seed = randdd.randInt();
    

    int const N(3000);   //number of vertices
    int Niter = 1000;
    double t(0);  
    
    int Exit = 0; 
    double r1; 
   
    
    double par = 1;
    
    double eps = 0.001;
	double gam = 1.5;
	double epsPow1mGam = pow(eps,1-gam);
   
/*   
    vector <double> act(N);
    vector <double> cum_act(N);
    double sum_act = 0;
    vector <double> cum_attrac(N);
    vector <double> attrac(N);
    double sum_attrac = 0;
   
    for (int i = 0; i < N; i++)
    {
        r1 = randdd.randDblExc();
        
        act[i] = pow(r1*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
        attrac[i] = 4;
        sum_act += act[i];
        cum_act[i] = sum_act;
        sum_attrac += attrac[i];
        cum_attrac[i] = sum_attrac;        
    }    
*/


    #pragma omp parallel for reduction(+:Exit) reduction(+:t) schedule(dynamic)
    for (int iter = 1; iter <= Niter; iter++)
    {
        MTRand mtrand1(iter+Seed);
	    vector <int> opinion(N);
        double r;
        int sum = 0;
        
        vector <double> act(N);
		vector <double> cum_act(N);
		double sum_act = 0;
		vector <double> cum_attrac(N);
		vector <double> attrac(N);
		double sum_attrac = 0;
	   
		for (int i = 0; i < N; i++)
		{
		    r = mtrand1.randDblExc();
		    
		    act[i] = pow(r*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
		    attrac[i] = 4;
		    sum_act += act[i];
		    cum_act[i] = sum_act;
		    sum_attrac += attrac[i];
		    cum_attrac[i] = sum_attrac;        
		}
        
        
        for (int i = 0; i < N/2; i++)
        {              
            opinion[i] = 1;
            sum++;                    
        }                   

        while (sum != 0 && sum != N)
        {
            r = mtrand1.randDblExc();
            int j = r*(N-1);
            
            r = r*sum_act;

			if (r > cum_act[j])
			{
				while (r > cum_act[j])
				{
					j++;								
				}
			}

			else
			{							
				while (r < cum_act[j-1] && j > 0)
				{
					j -= 1;
				}							
			}
            
            int Nguy = mtrand1.randInt(N-1);
            r = mtrand1.randDblExc();
            
            if (r < par)
            {
                sum += opinion[Nguy]-opinion[j];
                opinion[j] = opinion[Nguy];      // par = 1 voter model
            }
            else
            {
                sum += opinion[j]-opinion[Nguy];
                opinion[Nguy] = opinion[j];      // par = 0 invasion process
            }  

            r = mtrand1.randDblExc();               
            t += -log(1-r)/sum_act; 
            //t += 1./(N1*mu1);                                           
        }

        if (sum == N)
        {
            Exit++;
        }
    }
    
    cout << Exit/(double)Niter << " " << t/Niter << " Niter = " << Niter << " N = " << N << endl;
  
      
}


