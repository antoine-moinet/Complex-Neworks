#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include "MersenneTwister.h"
//#include <boost/math/special_functions/erf.hpp>
using namespace std;



int main()
{   
   
        
    MTRand mtrand1;

    int const N(100000);
    int const Niter(100);
    double t(1540);
    double tau;
    double cmax = 1;
    
    double alpha(0.9);
    unsigned long long int dmax(0);
    unsigned long long int ki = 0;
    double r(0.5);
    double r_avg(0);
    
    double t1 = 6.;
    double t0 = 0;
    
    double deltaT(t1/100.);
    int const N1 = (t1-t0)/deltaT+1;
    
    double Dmax(0);
    double eps(0.001);
    double betha(0.1);
   
    
    double kiMax = 0;
    double t_perc = 0;

    double *act = new double[N];
    int *count = new int[N];
    double *Time = new double[N];
    int *i_isat = new int[N];
    double *big2 = new double[N1];
    double *big = new double[N1];
    double *ki_s = new double[N1];
    
    vector< vector<int> > cluster(N);
   
    double A = betha/(pow(eps,-betha)-pow(cmax,-betha));
    
    
    
    double r2(0);
    
    for (int v = 0; v < N1; v++)
    {
        big[v] = 0;
        big2[v] = 0;
        ki_s[v] = 0;
    }

    for (int i = 0; i < N; i++)
    {
        r = mtrand1.randDblExc();
        //act[i] = pow(1-r,-1/betha);
        //act[i] = pow(1-r*betha/A,-1/betha);
        act[i] = pow(pow(eps,-betha)-r*betha/A,-1/betha);
        //act[i] = 1;
        r = mtrand1.randDblExc();
        //Time[i] = (pow(1-r,-1/alpha)-1)/act[i];
        Time[i] = -log(1-r)/act[i];
        i_isat[i] = i;
        cluster[i].push_back(i);
    }
    
    
    for (int iter = 1; iter <= Niter; iter++)
    {
    	for (int i = 0; i < N; i++)
		{
		    cluster[i].clear();
		    r = mtrand1.randDblExc();
		    //Time[i] = (pow(1-r,-1/alpha)-1)/act[i];
		    Time[i] = -log(1-r)/act[i];
		    i_isat[i] = i;
		    cluster[i].push_back(i);
		}
    	
    	for (int v = 0; v < N1; v++)
		{	    
		    //cout << v << endl;
		    
		    	
		    for (int compteur = 0 ; compteur < N ; compteur++)
		    {   
		        
		        while (Time[compteur] < v*deltaT+t0)
			    {		          
			        int n = mtrand1.randInt(N-1);
			        while (n == compteur)
			        {
			            n = mtrand1.randInt(N-1);
			        }
			        if (i_isat[compteur] != i_isat[n])
			        {
			        	int ici = i_isat[n];
			        	int la = i_isat[compteur];
						if (cluster[ici] > cluster[la])
						{
							ici = i_isat[compteur];
							la = i_isat[n];
						}
			        	
			        	for (int l = 1; l <= cluster[ici].size(); l++)
						{		            		
							cluster[la].push_back(cluster[ici][cluster[ici].size()-l]);
							i_isat[cluster[ici][cluster[ici].size()-l]] = la;				
						}

			        	//cluster[la].insert(cluster[la].end(), cluster[ici].begin(), cluster[ici].end());
			        	cluster[ici].clear();
			        				            	
			        }
			        r = mtrand1.randDblExc();
		        	//Time[compteur] += (pow(1-r,-1/alpha)-1)/act[compteur];
		        	Time[compteur] += -log(1-r)/act[compteur];
		        	
			    }
		    }
		    
		   
		    ki = 0;
		    dmax = 0; 		   
/*	  
		    for (int k = 0; k < N; k++)
			{   			   	
			    if (cluster[k].size() > dmax)
			    {
			        dmax = cluster[k].size();			        			        
			    }	
			    ki += cluster[k].size()*(cluster[k].size()-1);    
			}	
			
			ki -= dmax*(dmax-1);	
*/			
			int number_of_dmax = 0;
			for (int k = 0; k < N; k++)
			{   			   	
			    if (cluster[k].size() > dmax)
			    {
			        dmax = cluster[k].size();			        			        
			    }			    
		        ki += cluster[k].size()*cluster[k].size(); 				        			          		        			             				       
			}
			for (int k = 0; k < N; k++)
			{   			   	
			    if (cluster[k].size() == dmax)
			    {
			        number_of_dmax++;			        			        
			    } 			        			             				       
			}			
			ki -= dmax*dmax*number_of_dmax;			
					
			ki_s[v] += ki/((double)Niter*N);
		   
		    big[v] += dmax/(double)Niter;
		    big2[v] += dmax*dmax/(double)Niter;
		}
	}
        
    for (int v = 0; v < N1; v++)
	{
		/*if (ki_s[v] > kiMax)
        {
        	kiMax = ki_s[v];
        	t_perc = v*deltaT+t0;
        }*/
        if (big2[v]-pow(big[v],2.) > kiMax)
        {
        	kiMax = big2[v]-pow(big[v],2.);
        	t_perc = v*deltaT+t0;
        }
	}
        
        

    ofstream files("BigClustFast107.dat");
    for (int v = 1; v < N1; v++)
	{
		files << v*deltaT+t0 << " " << big[v]/(double)N << " " << pow((big2[v]-pow(big[v],2.))/((double)N*N) , 0.5) << " " << ki_s[v] << endl;
		//files << v*deltaT+t0 << " " << big[v]/(double)N << " " << ki_s[v] << endl;
	}
        
    files.close();   

       
       
    
    
    cout << t_perc << endl;
    
    delete [] count;
    delete [] act; 
    delete [] Time;
    delete [] i_isat;  
    delete [] big2;  
    delete [] big; 
    delete [] ki_s;
}



