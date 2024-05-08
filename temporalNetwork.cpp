#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "MersenneTwister.h"
#include <vector>
#include <omp.h>
using namespace std;


int main()
{       
   
    MTRand mtrand1;

    int NumPoints,N,r_avg2,r_avg;                        
    double act,alpha,t,alpha2;


    cout << "number of nodes?" << endl;
    cin >> N;
    cin.ignore();
    cout << "time?" << endl;
    cin >> t;
    cin.ignore();
    cout << "number of points?" << endl;
    cin >> NumPoints;
    cin.ignore();
    cout << "alpha?" << endl;
    cin >> alpha;
    cin.ignore();
    cout << "alpha2?" << endl;
    cin >> alpha2;
    cin.ignore();
    cout << "activity?" << endl;
    cin >> act;
    cin.ignore();
    
    cout << endl;
    

    vector<int> avg_r(NumPoints,0);
    vector<int> avg_r2(NumPoints,0);

/*
    A = betha/(pow(eps,-betha)-1);       

    for (int i = 0; i < N; i++)    ///////////////// assign a power law distributed activity c between epsilon and 1
    {
        r = mtrand1.randDblExc();            
        act[i] = pow(pow(eps,-betha)-r*betha/A,-1/betha);            
    }        
*/   
        
    for (int i = 0 ; i < avg_r.size() ; i++)
    {
      
        r_avg = 0;
        r_avg2 = 0;
        double tmax = (i+1)*t/NumPoints;
        int seed = mtrand1.randInt();
        cout << i << endl;

        #pragma omp parallel for reduction(+:r_avg) reduction(+:r_avg2) schedule(dynamic)
        for (int j = 0; j < N; j++)
        {
            double t1,t2,t3,r;
            MTRand mtrand2(j+seed);

            r = mtrand2.randDblExc();
            t1 = (pow(1-r,-1/alpha)-1)/act;
            t2 = 0;
            r = mtrand2.randDblExc();
            //t3 = (pow(1-r,-1/alpha)-1)/act;
			t3 = (exp(   pow(   pow(log(2),alpha2)  /(1-r),1/alpha2)   )-2)/act;
			
            while (t1 < tmax)
            {
                r_avg++;
                while (t2 < t1)
                {
                    r = mtrand2.randDblExc();
                    t2 += (pow(1-r,-1/alpha)-1)/act;
                }
                t1 = t2;
                t2 = 0;
            }

            while (t3 < tmax)
            {
                r_avg2++;
                r = mtrand2.randDblExc();
                //t3 += (pow(1-r,-1/alpha)-1)/act;
                t3 += (exp(pow(pow(log(2),alpha2)/(1-r),1/alpha2))-2)/act;
            }

        }   

        avg_r[i] = r_avg;
        avg_r2[i] = r_avg2;
    }
    
  
    ofstream opa ("RaverageNonRenewal.dat");

    for (int k = 0 ; k < avg_r.size(); k++)
    {
        opa << (k+1)*t/NumPoints << " " << avg_r[k]/(double)N << " " << avg_r2[k]/(double)N << endl;
        //opa << k << " " << pp[k] << " " << exp(-r_avg+k*log(r_avg)-logsum(k))  << endl;
        //opa << k << " " << p1[k] << " " << p1[k+(int)Dravg] - pow(ta,alpha)/gam*(p[k+(int)Dravg-1]-p[k+(int)Dravg]) << endl;
    }
    opa.close();    
   
}



