#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "MersenneTwister.h"
#include <vector>
#include <omp.h>
//#include <boost/math/special_functions/erf.hpp>

using namespace std;

double logsum(int n)
{
    if (n == 0)
    {
        return 0;
    }
    return log(n)+logsum(n-1);
}

int main()
{       
   
    MTRand mtrand1;

    int const N(1000000);
    int const Niter = 50;
    double t(100000);
    //double tau;
    double ta(1000);
    
    double alpha(0.5);
    double dmax(0);
    double r(0.5);
    double r_avg(0);
    double deltaT(1);
    double sum(0);
    double Dmax(0);
    double eps(0.001);
    double betha(3);
    double gamma = pow(tgamma(1-alpha),1/alpha);
    double A;
    double c_0 = 1;
    double delta;
    int sign;
    int newSign;
    double t_perc;

    double *act = new double[N];
    int *count = new int[N];

    double r2(0);   
    int nonZero = 0;
    double tmax = 10000;
    
    //A = betha/(pow(eps,-betha)-1);     
    
    for (int i = 0; i < N; i++)                    ////////////////////////  aged degree distribution
    {
        r = mtrand1.randDblExc();
        //act[i] = c_0*pow(1-r,-1/betha);
        //act[i] = pow(pow(eps,-betha)-r*betha/A,-1/betha);
        //act[i] = pow(1-r*(0.1*u+0.01)/A,-1/(0.1*u+0.01));
        act[i] = 0;     
        //act[i] = tmax*pow(r,1/(1+betha));  
    }        

    #pragma omp parallel for schedule(dynamic)
    for (int compteur = 0 ; compteur < N ; compteur++)
	{
    	MTRand mtrand2(compteur);
    	double tau;
    	r = mtrand2.randDblExc();
    	count[compteur] = 0;
    	//tau = -log(1-r)/act[compteur];
    	tau = pow(1-r,-1/alpha)-1;
    	//tau = pow(boost::math::erf_inv(1-r),-2)/act[compteur];
    	
   	 	while (tau < t + act[compteur])
    	{
        		if (tau > act[compteur])
        		{
            		count[compteur] += 1;
        		}                
        		r = mtrand2.randDblExc();
        		//tau += pow(boost::math::erf_inv(1-r),-2)/act[compteur];
        		tau += pow(1-r,-1/alpha)-1;
        		//tau += -log(1-r)/act[compteur];

    	}
/*    	
    	if (count[compteur] > 0)
    	{
    		nonZero += 1;
    	}

    	r_avg += count[compteur]/(double)N;
*/	
	}  
	//cout << r_avg << " " << nonZero << endl;      	  
		        
    int *G = new int[N];

    for (int k = 0; k < N; k++)
    {
        G[k] = 0;
    }
    
	#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < N; k++)
    {
        MTRand mtrand3(k);
        int j = count[k];
        while (j != 0)
        {
            int n = mtrand3.randInt(N-1);
            while (n == k)
            {
                n = mtrand3.randInt(N-1);
            }

            G[k] += 1;
            G[n] += 1;                
            j -= 1;                               
        }
    }

    

    for (int k = 0; k < N; k++)
    {
        if (G[k] > Dmax)
        {
            Dmax = G[k];
        }
    }

    int const D(Dmax+1);
    vector<double>  p(D);


    for (int i = 0 ; i < D ; i++)
    {
        p[i] = 0;
    }
	
	
	#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; i++)
    {
        int j = 0;
        while (G[i] != j)
        {
            j++;
        }
        p[j] += 1/(double)N;
    }       


    ofstream op ("AgeDegDistr.dat");

    for (int compteur = 0 ; compteur < D; compteur++)
    {
        op << compteur << " " << p[compteur] << endl;
    }
    op.close();
    
    p.clear();

/*
    double r_avg2(0);
    
    for (int compteur = 0 ; compteur < N ; compteur++)              //////////////////// non-aged degree distribution
	{
	    r = mtrand1.randDblExc();
	    count[compteur] = 0;
	    //tau = -log(1-r)/act[compteur];
	    tau = (pow(1-r,-1/alpha)-1)/(act[compteur]*gamma);
	    while (tau < t)
	    {            
	        count[compteur] += 1;                
	        r = mtrand1.randDblExc();
	        tau += (pow(1-r,-1/alpha)-1)/(act[compteur]*gamma);
	        //tau += -log(1-r)/act[compteur];
	    }
	    r_avg2 += count[compteur]/(double)N;		    
	}
    
    for (int k = 0; k < N; k++)
	{
	    G[k] = 0;
	}	

	for (int k = 0; k < N; k++)
	{
	    int j = count[k];
	    while (j != 0)
	    {
	        int n = mtrand1.randInt(N-1);
	        while (n == k)
	        {
	            n = mtrand1.randInt(N-1);
	        }

	        G[k] += 1;
	        G[n] += 1;                
	        j -= 1;                               
	    }
	}

	double Dmax1 = 0;

	for (int k = 0; k < N; k++)
	{
	    if (G[k] > Dmax1)
	    {
	        Dmax1 = G[k];
	    }
	}	
	
    vector<double> p1(Dmax1+1);			

	for (int i = 0; i < N; i++)
	{
	    int j = 0;
	    while (G[i] != j)
	    {
	        j++;
	    }
	    p1[j] += 1/(double)N;
	}    
    
    betha = betha-alpha;
    
    for (int i = 0; i < N; i++)                        ///////////////////// non-aged degree distribution with betha --> betha-alpha
    {
        r = mtrand1.randDblExc();
        act[i] = c_0*pow(1-r,-1/betha);
        //act[i] = pow(pow(eps,-betha)-r*betha/A,-1/betha);
        //act[i] = pow(1-r*(0.1*u+0.01)/A,-1/(0.1*u+0.01));
        //act[i] = c_0;       
    }
    
    double r_avg22 = 0;
    
    for (int l = 1; l <= Niter; l++)
    {
    	cout << l << endl;
    	
    	for (int compteur = 0 ; compteur < N ; compteur++)
		{
		    r = mtrand1.randDblExc();
		    count[compteur] = 0;
		    //tau = -log(1-r)/act[compteur];
		    tau = (pow(1-r,-1/alpha)-1)/(act[compteur]*gamma);
		   		                
		    while (tau < t)
		    {
		    	count[compteur] += 1;		    	              
				r = mtrand1.randDblExc();
				tau += (pow(1-r,-1/alpha)-1)/(act[compteur]*gamma);
				//tau += -log(1-r)/act[compteur];
		    }		        
		    r_avg22 += count[compteur]/((double)N*Niter);	    
		}       
        
		for (int k = 0; k < N; k++)
		{
		    G[k] = 0;
		}
		

		for (int k = 0; k < N; k++)
		{
		    int j = count[k];
		    while (j != 0)
		    {
		        int n = mtrand1.randInt(N-1);
		        while (n == k)
		        {
		            n = mtrand1.randInt(N-1);
		        }

		        G[k] += 1;
		        G[n] += 1;                
		        j -= 1;                               
		    }
		}

		Dmax = 0;

		for (int k = 0; k < N; k++)
		{
		    if (G[k] > Dmax)
		    {
		        Dmax = G[k];
		    }
		}	
		
	    while (p.size() < Dmax+1)
	    {
	    	p.push_back(0);
	    }			

		for (int i = 0; i < N; i++)
		{
		    int j = 0;
		    while (G[i] != j)
		    {
		        j++;
		    }
		    p[j] += 1/((double)N*Niter);
		}
    }

    double Dravg = r_avg2-r_avg;
    double Dravg2 = r_avg22-r_avg;

    double gam = tgamma(1+alpha);    //gamma(1+alpha)

    ofstream opa ("degDistr2.dat");
    
    double kmin = p.size();
    if (Dmax1 < kmin)
    {
    	kmin = Dmax1;
    }

    for (int k = 0 ; k < kmin-(int)Dravg; k++)
    {
        opa << k << " " << p1[k] << " " << p1[k+(int)Dravg] - (betha+alpha)/betha*pow(ta,alpha)/gam*(p[k+(int)Dravg2-1]-p[k+(int)Dravg2]) << " " << p[k] << endl;
        //opa << k << " " << pp[k] << " " << exp(-r_avg+k*log(r_avg)-logsum(k))  << endl;
        //opa << k << " " << p1[k] << " " << p1[k+(int)Dravg] - pow(ta,alpha)/gam*(p[k+(int)Dravg-1]-p[k+(int)Dravg]) << endl;
    }
    opa.close();
    
 */  
   
    delete [] G;       
    delete [] count;    
    delete [] act; 
    //cout << gam << " " << gamma << endl;  
}



