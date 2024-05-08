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

    int const N(30);
    int const Niter(100000);   //number of vertices
    double t(100000);  
    double deltaT(1000);  
    double alpha(0.3);
    int const N1((int)(t/deltaT));
    //int N1 = 100;
    
  
    int u;
    double tau;
    int Nguy;
    double r;
    int rho(0);
   
    double *Time = new double[N];
    int *opinion = new int[N];
    int *rhoo = new int[N1+1];             //////////////////////// density as a function of time      
    double *count = new double[N1+1];    ////////////////////////// density as a function of time     
    int rhoi;


    
                                      //////////////////////////// density as a function of time 
    for (int k = 0; k < N1+1; k++)
    {        
        count[k] = Niter;
        rhoo[k] = 0;
    }

    for (int m = 0; m < Niter; m++)
    {
        //cout << m << endl;
        rho = 0;
        while (rho == 0 || rho == N)
        {
            rho = 0;

            for (int k = 0; k < N; k++)
            {        
                opinion[k] = 0;
                r = mtrand1.rand();
                if (r < 0.5)  /////////////////////////////////////////////// random initial conditions
                //if (k < N/2)  ///////////////////////////////////////////////// rho_0 = 0.5
                {
                    opinion[k] = 1;
                }
            }

            for (int k = 0; k < N; k++)
            {
                rho += opinion[k];
            }
        }

        rhoi = rho*(N-rho);
        rhoo[0] += rhoi;
        
        for (int k = 0; k < N; k++)
        {        
            Time[k] = 0;
            r = mtrand1.rand();
            Time[k] += pow(1-r,-1/alpha)-1;  
            //Time[k] = -log(1-r);      
        }

        tau = Time[0];
        u = 0;

        for (int k = 0; k < N; k++)
        {                
            if (Time[k] < tau)
            {
                tau = Time[k];
                u = k;
            }          
        }
        
        for (int i = 1; i <= N1; i++)
        {            
            while (tau < i*deltaT)
            {
                r = mtrand1.rand();
                Time[u] += pow(1-r,-1/alpha)-1;
                //Time[u] += -log(1-r);
                tau = Time[u];

                Nguy = mtrand1.randInt(N-1);
                while (Nguy == u)
                {
                    Nguy = mtrand1.randInt(N-1);                    
                } 
                                   
               
                if (opinion[u] == 0)/////////////////////////////////////// voter density of interfaces
                {
                    rhoi += (opinion[Nguy]-opinion[u])*(N-2*rho-1);
                }

                else
                {
                    rhoi += (opinion[u]-opinion[Nguy])*(2*rho-N-1);
                }

                rho += opinion[Nguy]-opinion[u];/////////////////////////// voter 
                opinion[u] = opinion[Nguy];                
           
/*              
                if (opinion[Nguy] == 0)/////////////////////////////////////////////// moran density of interfaces
                {
                    rhoi += (opinion[u]-opinion[Nguy])*(N-2*rho-1);
                }

                else
                {
                    rhoi += (opinion[Nguy]-opinion[u])*(2*rho-N-1);
                }

                rho += opinion[u]-opinion[Nguy];////////////////////////// moran regular density
                opinion[Nguy] = opinion[u]; 
*/

                for (int k = 0; k < N; k++)
                {                
                    if (Time[k] < tau)
                    {
                        tau = Time[k];
                        u = k;
                    }          
                }           
            }



            if (rho == 0 || rho == N)
            {
                count[i] -= 1;
            }

            //if (rho != 0 && rho != N)
            //{
            rhoo[i] += rhoi;//////////////////// density of surfaces
                //rhoo[i] += rho;///////////////////// regular density
            //}  
                 
        }       
    }  
   
    ofstream files ("NonMarkovVoter3-1000.dat");
    
    for (int i = 0; i <= N1; i++)
    {
        //files << i*deltaT << " " << 2*rhoo[i]/(N*(N-1)*(count[i]+0.0001)) << endl;//////////// density of surfaces
        files << i*deltaT << " " << 2*rhoo[i]/(N*(N-1)*(double)Niter) << endl;//////////// density of surfaces
        //files << i*deltaT << " " << rhoo[i]/(N*count[i]) << endl;//////////////////// regular density
        cout << count[i] << " " << i << endl;
    } 

    files.close();

    delete[] opinion;
    delete[] Time;   
    //delete[] rhoo;   
    //delete[] count;
   
    opinion = 0;
    Time = 0;
   
}


