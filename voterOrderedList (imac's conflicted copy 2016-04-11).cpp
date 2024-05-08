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
   

    int const N(30);   //number of vertices
    int const Niter(1000);
    double t(100000);  
    double deltaT(100000);  
    double alpha(0.5);
  
    int j;
    int i(1);
    int u;
    int Nguy;
    double r;
    int rho(0);
    double *Time = new double[N];
    int *opinion = new int[N];
    double *consensusTime = new double[Niter];
    vector<double> list;
    vector<int> listIndex;
    vector<double> vide;
    vector<int> videIndex;
    double a;

    
    for (int h = 0; h < Niter; h++)
    {
        rho = 0;
        i = 1;

        for (int k = 0; k < N; k++)
        {
            Time[k] = 0;
            opinion[k] = 0;
            r = mtrand1.rand();
            if (r < 0.5)
            {
                opinion[k] = 1;
            }
        }

        for (int k = 0; k < N; k++)
        {
            rho += opinion[k];
        }

        //ofstream files ("VoterOneRun.dat");
        //files << 0 << " " << rho/(double)N << endl;
        //files.close();

        for (int k = 0; k < N; k++)
        {        
            r = mtrand1.rand();
            Time[k] += pow(1-r,-1/alpha)-1;        
        }
        
        while (rho != 0 && rho != N)
        {
            list = vide;
            listIndex = videIndex;
            for (int k = 0; k < N; k++)
            {
                while (Time[k] < i*deltaT)
                {
                    list.push_back(Time[k]);
                    listIndex.push_back(k);
                    r = mtrand1.rand();
                    Time[k] += pow(1-r,-1/alpha)-1;
                }            
            }

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

            for (int k = 0; k < list.size() ; k++)
            {   
                Nguy = mtrand1.randInt(N-1);

                rho += opinion[Nguy]-opinion[listIndex[k]];///////////////////////// voter                
                opinion[listIndex[k]] = opinion[Nguy];

                //rho += opinion[listIndex[k]]-opinion[Nguy];///////////////////////// moran
                //opinion[Nguy] = opinion[listIndex[k]];  
            }

            //ofstream files ("VoterOneRun.dat",ios::app);
            //files << i*deltaT << " " << rho/(double)N << endl;
            //files.close();
            i++;
        }

        consensusTime[h] = i*deltaT;
    }

    double *p = new double[10];



    for (i = 1; i <= 10; i++)
    {
        p[i] = 0;
    }

    for (int h = 0; h < Niter; h++)
    {
        i = 1;
        while (consensusTime[h] < i*deltaT && i < 10)
        {
            i++;
        }
        p[i] += 1/(double)Niter;
    }
    
    ofstream files ("VoterOneRun.dat");
    for (i = 1; i <= 10; i++)
    {
        files << i*deltaT << " " << p[i] << endl;
    }
    
    files.close();


    delete[] opinion;
    delete[] Time;
    delete[] consensusTime; 
    delete[] p;    
}


