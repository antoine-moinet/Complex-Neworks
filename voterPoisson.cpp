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
    int const Niter(100000);   
    double t(100);  
    double deltaT(1);  
    double alpha(0.3);
    int const N1((int)(t/deltaT));
    
    
  
    int j = 1;
    int firstGuy;
    int secondGuy;
    double r;
    int rho(0);
    int rho_0(0);
    double tau = 0;
    int *opinion = new int[N];
    int *rhoo = new int[N1+1];
    int *surfRho = new int[N1+1];            
    int rhoi;

    for (int k = 0; k < N1; k++)
    {
        rhoo[k] = 0;
        surfRho[k] = 0;
    }
  

    for (int h = 0; h < Niter; h++)
    {
        //cout << h << endl;
        j = 1;
        rho_0 = 0;
        tau = 0;
        while (rho_0 == 0 || rho_0 == N)
        {
            rho_0 = 0;

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
                rho_0 += opinion[k];
            }
        }

        rho = rho_0;
        rhoi = rho*(N-rho);
        rhoo[0] += rho;
        surfRho[0] += rhoi;

        r = mtrand1.rand();

        tau += -log(1-r)/(double)N;
        firstGuy = mtrand1.randInt(N-1);
        secondGuy = mtrand1.randInt(N-1);

        while (firstGuy == secondGuy)
        {
            secondGuy = mtrand1.randInt(N-1);                    
        } 

        while (j < N1)
        {
            while (tau < j*deltaT)
            {
/*               
                if (opinion[firstGuy] == 0)/////////////////////////////////////// voter density of interfaces
                {
                    rhoi += (opinion[secondGuy]-opinion[firstGuy])*(N-2*rho-1);
                }

                else
                {
                    rhoi += (opinion[firstGuy]-opinion[secondGuy])*(2*rho-N-1);
                }

                rho += opinion[secondGuy]-opinion[firstGuy];/////////////////////////// voter 
                opinion[firstGuy] = opinion[secondGuy];                
*/
                
                if (opinion[secondGuy] == 0)/////////////////////////////////////////////// moran density of interfaces
                {
                    rhoi += (opinion[firstGuy]-opinion[secondGuy])*(N-2*rho-1);
                }

                else
                {
                    rhoi += (opinion[secondGuy]-opinion[firstGuy])*(2*rho-N-1);
                }

                rho += opinion[firstGuy]-opinion[secondGuy];////////////////////////// moran regular density
                opinion[secondGuy] = opinion[firstGuy]; 
              

                r = mtrand1.rand();
                tau += -log(1-r)/(double)N;
                firstGuy = mtrand1.randInt(N-1);
                secondGuy = mtrand1.randInt(N-1);

                while (firstGuy == secondGuy)
                {
                    secondGuy = mtrand1.randInt(N-1);                    
                } 
                
            }

            surfRho[j] += rhoi;
           
            rhoo[j] += rho;
            j++;
        }
        
    }

    ofstream files ("voterPoisson.dat");
    
    for (int i = 0; i < N1; i++)
    {        
        files << i*deltaT << " " << 2*surfRho[i]/(N*(N-1)*(double)Niter) << " " << surfRho[i] << endl;
        //files << i*deltaT << " " << rhoo[i]/(N*count[i]) << endl;//////////////////// regular density
        //cout << count[i] << " " << i << endl;
    } 

    files.close();

    delete[] opinion;
    delete[] surfRho;
    delete[] rhoo;    
}




