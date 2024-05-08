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
    
  
    int j;
    int u;
    int Nguy;
    double r;
    int rho(0);
    int rho_0(0);
    double *Time = new double[N];
    int *opinion = new int[N];
    int *opinion_0 = new int[N];
    int *rhoo = new int[N1+1];             //////////////////////// density as a function of time      
    double *count = new double[N1+1];    ////////////////////////// density as a function of time  
    //vector<int> rhoo(N1+1,0);                     /////////////////////////// density as a function of the number of updates
    //vector<double> count(N1+1,Niter);                    /////////////////////////// density as a function of the number of updates
    vector<double> list;
    vector<int> listIndex;
    vector<double> vide;
    vector<int> videIndex;
    double a;
    int rhoi;
    int taille(N*pow(t,alpha));

    
                                      //////////////////////////// density as a function of time 
    for (int k = 0; k < N1+1; k++)
    {        
        count[k] = Niter;
        rhoo[k] = 0;
    }

 /*   
    for (int k = 0; k < N; k++)
    {        
        opinion_0[k] = 0;
        //r = mtrand1.rand();
        //if (r < 0.5)  /////////////////////////////////////////////// random initial conditions
        if (k < N/2)  ///////////////////////////////////////////////// rho_0 = 0.5
        {
            opinion_0[k] = 1;
        }
    }

    for (int k = 0; k < N; k++)
    {
        rho_0 += opinion_0[k];
    }
*/
    //rhoo[0] = Niter*rho_0;/////////////////////////////////// regular density
    rhoo[0] = Niter*rho_0*(N-rho_0);/////////////////////////// density of interfaces
    count[0] = Niter;

    for (int m = 0; m < Niter; m++)
    {
        //cout << m << endl;
        rho_0 = 0;
        while (rho_0 == 0 || rho_0 == N)
        {
            rho_0 = 0;

            for (int k = 0; k < N; k++)
            {        
                opinion_0[k] = 0;
                r = mtrand1.rand();
                if (r < 0.5)  /////////////////////////////////////////////// random initial conditions
                //if (k < N/2)  ///////////////////////////////////////////////// rho_0 = 0.5
                {
                    opinion_0[k] = 1;
                }
            }

            for (int k = 0; k < N; k++)
            {
                rho_0 += opinion_0[k];
            }
        }

        rho = rho_0;
        rhoi = rho*(N-rho);
        rhoo[0] += rhoi;
        
        for (int k = 0; k < N; k++)
        {        
            opinion[k] = opinion_0[k];
            Time[k] = 0;
            r = mtrand1.rand();
            Time[k] += pow(1-r,-1/alpha)-1;  
            //Time[k] = -log(1-r);      
        }
        
        for (int i = 1; i <= N1; i++)
        {
            list = vide;
            listIndex = videIndex;
            //u = 1;

            //while (list.size() < N1)
            //{
            for (int k = 0; k < N; k++)
            {
                //while (Time[k] < u*deltaT)
                while (Time[k] < i*deltaT)
                {
                    list.push_back(Time[k]);
                    listIndex.push_back(k);
                    r = mtrand1.rand();
                    Time[k] += pow(1-r,-1/alpha)-1;
                    //Time[k] = -log(1-r); 
                }            
            }
                //u++;

            //}
            //cout << "1ere boucle while passee avec t = " << u*deltaT <<  endl;

            for (int k = 1; k < list.size() ; k++)
            {   
                u = k;
                while (u > 0 && list[u] < list[u-1])
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

            //for (int k = 0; k < N1 ; k++)    ///////////////////// density as a function of the number of updates
            for (int k = 0; k < list.size() ; k++)
            {                  

                Nguy = mtrand1.randInt(N-1);
                while (Nguy == listIndex[k])
                {
                    Nguy = mtrand1.randInt(N-1);                    
                } 
                               
/*
                if (opinion[listIndex[k]] == 0)/////////////////////////////////////// voter density of interfaces
                {
                    rhoi += (opinion[Nguy]-opinion[listIndex[k]])*(N-2*rho-1);
                }

                else
                {
                    rhoi += (opinion[listIndex[k]]-opinion[Nguy])*(2*rho-N-1);
                }

                rho += opinion[Nguy]-opinion[listIndex[k]];/////////////////////////// voter 
                opinion[listIndex[k]] = opinion[Nguy];                
*/
       
                if (opinion[Nguy] == 0)/////////////////////////////////////////////// moran density of interfaces
                {
                    rhoi += (opinion[listIndex[k]]-opinion[Nguy])*(N-2*rho-1);
                }

                else
                {
                    rhoi += (opinion[Nguy]-opinion[listIndex[k]])*(2*rho-N-1);
                }

                rho += opinion[listIndex[k]]-opinion[Nguy];////////////////////////// moran regular density
                opinion[Nguy] = opinion[listIndex[k]]; 

/*                
                if (rho == 0 || rho == N)                          ///////////////////// density as a function of the number of updates
                {
                    count[k+1] -= 1;
                }

                if (rho != 0 && rho != N)
                {
                    rhoo[k+1] += rhoi;//////////////////// density of surfaces
                    //rhoo[i] += rho;///////////////////// regular density
                }                
*/                                        
            }  
       
                                                                ////////////////////////// density as a function of time 
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
    delete[] opinion_0;
    opinion = 0;
    Time = 0;
    opinion_0 = 0;
}


