#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "MersenneTwister.h"
using namespace std;



int main()
{
    MTRand mtrand1;

    int const N(100000);
    double ta(1000);
    double t(1000000);
    double tau;
    double alpha(0.5);
    double r(0.5);
    double deltaT(0.01);
    int const n(t/deltaT);
    double c(10);
    double tmax(0);
    double sinc = pow(sin(4*atan(1)*alpha)/atan(1)/4/alpha,1/alpha);
    double c2 = 1/sinc/(ta+1/c/sinc);
    cout << c2 << endl;
    
    double theta2 = 1./(alpha/(1-alpha)*pow(c*ta,alpha-1)*(pow(c*ta+1,1-alpha)-1)+1);    
    
    double *h = new double[n];
    double *h2 = new double[n];
    
    for (int k = 1; k <= n; k++)
    {
        h[k-1] = 0;
        h2[k-1] = 0;
    }

    for (int compteur = 0 ; compteur < N ; compteur++)
    {
        tau = 0.0;

        while (tau <= ta)
        {
            r = mtrand1.randDblExc();
            tau += (pow(1-r,-1/alpha)-1)/c;            
        }

        double a = (tau-ta)/deltaT;
        if (tau-ta > tmax)
        {
            tmax = tau-ta;
        }
		
		if (tau-ta)
		{
			
		}

        if (tau-ta < t)
        {
            h[(int)a] += 1./((double)N*deltaT);
        }    
        
        r = mtrand1.randDblExc();
        if (r < theta2)
        {
        	r = mtrand1.randDblExc();
        	tau = ta*pow(1-r,-1/alpha);  
        }
        else
        {
        	r = mtrand1.randDblExc();
        	tau = (pow((pow(c*ta+1,1-alpha)-1)*r+1,1/(1-alpha))-1)/c;
        }
        //tau = (pow(1-r,-1/alpha)-1)/c2;            
        
        if (tau < t)
        {
            h2[(int)(tau/deltaT)] += 1./((double)N*deltaT);
        }                
    }

    double sum = 0;
    for (int k = 1; k <= n; k++)
    {
        sum += h[k-1]*deltaT;
    }
    cout << sum << endl;
    cout << tmax << endl;
    
    int Numpoints = 100;
  
    double a = pow(deltaT/t,1./(1-Numpoints));
    
    
    double *g = new double[Numpoints];
    double *g2 = new double[Numpoints];

    for (int k = 0; k < Numpoints; k++)
    {
        g[k] = 0;
        g2[k] = 0;
        int count = 0;
        for (int i = (int)(t*pow(a,k-Numpoints)/deltaT); i < (int)(t*pow(a,k-Numpoints+1)/deltaT); i++)
        {
            count ++;
            g[k] += h[i];
        }
        g[k] = g[k]/(double)count;
        
        count = 0;
        for (int i = (int)(t*pow(a,k-Numpoints)/deltaT); i < (int)(t*pow(a,k-Numpoints+1)/deltaT); i++)
        {
            count ++;
            g2[k] += h2[i];
        }
        g2[k] = g2[k]/(double)count;
        //cout << g[k] << " " << k << endl;
    }

    ofstream files ("forwardT2.dat");
    for (int k = 1; k < 100; k++)
    {
        files << k*deltaT << " " << h[k] << " " << h2[k] << endl;
    }
  
    for (int k = (int)(Numpoints-1+log(100*deltaT/t)/log(a)); k < Numpoints; k++)
    {
        files << t*pow(a,k-Numpoints+1) << " " << g[k] << " " << g2[k] << endl;
    }

    files.close();
  
    delete[] h; 
    delete[] g;
   
}
