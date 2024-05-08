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

    int const N1(10000);   //number of vertices
    int t(0);    
    vector<double> pr(100,1/(double)N1);
    vector<double> prcum(101,1/(double)N1);
    int Nguy;
    double r;
    double gam(0.8);
    double eps(0.01);
    double som(0); 
    double sum(0);
    double tot(0);
    double tut(0);
    double par(0.2);  
    int const N(10000); 
    double *act = new double[N1]; 
    double *opinion = new double[N1];
    int *count = new int[N1];
    vector< vector<int> > Na(100);
    
    double mu1(0);
    double mu2(0);
    double mu_1(0); 
    double muOm(0);
    double Mupar(0);
    int cunt; 
    int j;

    double rho[100][1000];
    double rhoa[100];
    double Om[1000];
    double om1[1000];
    double om_1[1000];

    for (int k = 0; k < 100; k++)
    {
        tot += pow(eps*pow(1/eps,k/99.),-gam);
        tut += eps*pow(1/eps,k/99.);
        rhoa[k] = 0;
    }  
   
    for (int k = 0; k <= 100; k++)
    {   
        prcum[k] = som;
        som += pow(eps*pow(1/eps,k/99.),-gam)/tot;
        if (k < 100)
        {
            mu_1 += pow(eps*pow(1/eps,k/99.),-gam-1)/tot;
            mu1 += pow(eps*pow(1/eps,k/99.),-gam+1)/tot;
            mu2 += pow(eps*pow(1/eps,k/99.),-gam+2)/tot;
        }
    }  
    cout << mu_1 << " " << mu1 << " " << mu2 << endl;

    for (int k = 0; k < 100; k++)
    {           
        Mupar += pow(eps*pow(1/eps,k/99.),-gam)/(tot*(par*eps*pow(1/eps,k/99.)+(1-par)*mu1));        
    }  

    for (int k = 0; k < 100; k++)
    {   
        muOm += pow(eps*pow(1/eps,k/99.),-gam)*((1-par)*(mu1-eps*pow(1/eps,k/99.))+par*eps*pow(1/eps,k/99.))/(tot*Mupar*par*(par*eps*pow(1/eps,k/99.)+(1-par)*mu1));
    }  

    mu_1 = 0;
    mu1 = 0;
    mu2 = 0;

    for (int k = 0; k < 1000; k++)
    {   
        Om[k] = 0;
        om1[k] = 0;
        om_1[k] = 0;  
        for (int z = 0; z <= 99; z++)
        {
            rho[z][k] = 0;
        }       
    }  

    for (int i = 0; i < N1; i++)
    {
        r = mtrand1.rand();
        cunt = 0;
        while (prcum[cunt+1] < r)
        {
            cunt ++;   
        }
        
        Na[cunt].push_back(i); 
        count[i] = cunt;
        act[i] = eps*pow(1/eps,cunt/99.); 
        opinion[i] = 0;
    } 

    som = 0;   
     
    for (int k = 0; k <= 100; k++)
    {   
        prcum[k] = som;
        som += eps*pow(1/eps,k/99.)/tut;
    }  

    
    for (int i = 0; i < N1; i++)
    {
        mu_1 += 1./(N1*act[i]);
        mu1 += act[i]/(double)N1;
        mu2 += act[i]*act[i]/(double)N1;
    }
    
    cout << mu_1 << " " << mu1 << " " << mu2 << endl;

    /*for (int k = 0; k < 100; k++)
    {
        cout << Na[k].size() << " " << k << endl;
    }*/

    int p = 3;
    int a = 9;

    t = 0; 

    for (int i = 0; i < p*N1/21.; i++)
    {              
        opinion[i] = 1;
        rhoa[count[i]] += 1./Na[count[i]].size();         
    }

    for (int z = 0; z <= 99; z++)
    {        
        rho[z][0] = rhoa[z];       
    }   

    for (int k = 0; k < N; k++)
    {    
        t = 0;
        //cout << k << endl;

        for (int z = 0; z < 100; z++)
        {        
            rhoa[z] = 0;       
        }   

        for (int i = 0; i < N1; i++)
        {   
            opinion[i] = 0;

            if (i < p*N1/21.)
            {
                opinion[i] = 1;
                rhoa[count[i]] += 1./Na[count[i]].size();
            }    
        }

        while (t < a*N1)
        {
            r = mtrand1.randDblExc();
            cunt = 0;
            while (prcum[cunt+1] < r)
            {
                cunt ++;   
            }

            j = mtrand1.randInt(Na[cunt].size()-1);
            j = Na[cunt][j];
                       
            Nguy = mtrand1.randInt(N1-1);
            r = mtrand1.rand();
            if (r < par)
            {                              
                rhoa[cunt] += (opinion[Nguy]-opinion[j])/(double)Na[cunt].size();               
                opinion[j] = opinion[Nguy];      // par = 1 voter model
            }
            else
            {
                rhoa[count[Nguy]] += (opinion[j]-opinion[Nguy])/(double)Na[count[Nguy]].size();                        
                opinion[Nguy] = opinion[j];      // par = 0 invasion process
            }                              
                           
            t += 1; 
            if ((1000*t)%(a*N1) == 0 && 1000*t/(a*N1) < 1000)
            {                
                for (int z = 0; z <= 99; z++)
                {  
                    rho[z][1000*t/(a*N1)] += rhoa[z]/(double)N;
                }                
            }        
        }    
    } 

    for (int k = 0; k < 1000; k++)
    {
        for (int z = 0; z < 100; z++)
        {
            Om[k] += Na[z].size()/(double)N1*(par*muOm+(1-par)*act[Na[z][0]])*rho[z][k]/(par*act[Na[z][0]]+(1-par)*mu1);
            om1[k] += Na[z].size()/(double)N1*act[Na[z][0]]*rho[z][k]/mu1;
            om_1[k] += Na[z].size()/(double)N1*rho[z][k]/(act[Na[z][0]]*mu_1);
        }
    }

    ofstream fyless ("voter_rhoa.dat");
    for (int k = 0; k < 1000; k++)
    {
        for (int z = 0; z <= 9; z++)
        {  
            fyless << rho[11*z][k] << " ";
        } 
        fyless << Om[k] << " ";
        fyless << om1[k] << " ";
        fyless << om_1[k] << " "; 
        fyless << a*k/(1000*mu1) << endl;
    }
    fyless.close(); 

    cout << "wesh ma gueule" << endl;

    delete[] opinion;  
    delete[] act;     
}


