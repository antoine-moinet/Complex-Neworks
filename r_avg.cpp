#include <iostream>             // les 4 1eres lignes c'est pour appeler des librairies pour sauvegarder des fichiers et utiliser des fonctions mathematiques
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include "MersenneTwister.h"    // les deux c'est pour mettre un commentaire cette ligne te permet d'utiliser mersenne
//#include <boost/math/special_functions/erf.hpp>

using namespace std;

int main()           // ça c'est ton programme principal ça commence tout le temps comme ça
{   
    ofstream files ("r_avg.dat");
    files << 0 << " " << 0 << " " << 0 << endl;
    files.close();

    
    MTRand mtrand1;          // ça c'est pour pouvoir utiliser twister il faut le déclarer

    int const N(1000000);  // tu déclare toutes tes variables en leur assignant ou non une valeur les () ou = c'est pareil
    int N1 = 100;
    
    double deltaT = 0.3;
    double r_avg, tau, r, r2, r_i;
    double ta = 0;
    //double alpha = -2;
    double alpha = 0.9;
    double betha = 0.1;
    double eps = 0.001;
    double cmax = 1;
    
    double A = betha/(pow(eps,-betha)-pow(cmax,-betha));
    vector<int> count(N,0);
    vector<double> act(N);
    r_avg = 0;
    r2 = 0;
    r_i = 0;
    
    for (int j = 0; j < N; j++)         // une boucle for
    {
    	r = mtrand1.randDblExc();
        act[j] = pow(pow(eps,-betha)-r*betha/A,-1/betha);
        r = mtrand1.randDblExc();
        //tau = (pow(1-r,-1/alpha)-1)/act[j];
        tau = -log(1-r)/act[j];
        while (tau < 30)  // dans les conditions si tu veux mettre un égal il faut écrire ==, après tu as <= (inférieur ou égal) >= et != (différent de)
        {
            count[j]++;            
            r = mtrand1.randDblExc();
            //tau += pow(boost::math::erf_inv(1-r),-2); ///////////////// levy distr.
            //tau += pow(-log(r),-1/alpha);              //////////////////  levy-Weibull distr.
            //tau += (pow(1-r,-1/alpha)-1)/act[j];
            tau += -log(1-r)/act[j]; 
        }
        
        r_avg += count[j];
        r2 += count[j]*count[j];                
    }
    
    
    
    double avg_a = r_avg/N/30;  ///////// divided by t = 1
    double avg_a2 = (r2/N-r_avg/N)/pow(30,2.);
    
    double avg_a2_1 = r2/N/pow(30,2.);
    double avg_a2_2;
    
    r_avg = 0;
    r2 = 0;
    
    for (int j = 0; j < N; j++)         // une boucle for
    {
    	count[j] = 0;
    	for (int iter = 1; iter <= 100; iter++)
    	{
    		r = mtrand1.randDblExc();
		    //tau = (pow(1-r,-1/alpha)-1)/act[j];
		    tau = -log(1-r)/act[j];
		    while (tau < 30)  // dans les conditions si tu veux mettre un égal il faut écrire ==, après tu as <= (inférieur ou égal) >= et != (différent de)
		    {
		        count[j]++;            
		        r = mtrand1.randDblExc();
		        //tau += pow(boost::math::erf_inv(1-r),-2); ///////////////// levy distr.
		        //tau += pow(-log(r),-1/alpha);              //////////////////  levy-Weibull distr.
		        //tau += (pow(1-r,-1/alpha)-1)/act[j];
		        //tau += -log(1-r)*2; 
		        tau += -log(1-r)/act[j];
		    }	
    	}
    	
        
        r_avg = count[j]/100.;
        r2 += r_avg*r_avg;                
    }
    avg_a2_2 = r2/N/pow(30,2.);
    
 
    for (int i = 1; i <= N1; i++)         // une boucle for
    {
       
        r_avg = 0;
        r2 = 0;

        for (int j = 0; j < N; j++)         // une boucle for
        {
            
            
            r_i = 0;
            r = mtrand1.randDblExc();
            //tau = pow(boost::math::erf_inv(1-r),-2);
            //tau = pow(-log(r),-1/alpha);
            //tau = (pow(1-r,-1/alpha)-1)/act[j];  /////////// equilibrium process: first inter-event distributed with h -> alpha-1 
            tau = -log(1-r)/act[j];
            //tau = -log(1-r)*2;

            while (tau < ta + i*deltaT)  // dans les conditions si tu veux mettre un égal il faut écrire ==, après tu as <= (inférieur ou égal) >= et != (différent de)
            {
                if (tau > ta)
                {
                    r_i++;
                    
                    r = mtrand1.randDblExc();
                    //tau += pow(boost::math::erf_inv(1-r),-2); ///////////////// levy distr.
                    //tau += pow(-log(r),-1/alpha);              //////////////////  levy-Weibull distr.
                    //tau += (pow(1-r,-1/alpha)-1)/act[j];
                    tau += -log(1-r)/act[j];
                    //tau += -log(1-r)*2;                             
                }
            }
            r_avg += r_i;
            r2 += r_i*r_i;
        }

        ofstream fyl ("r_avg.dat",ios::app);
        fyl << i*deltaT << " " << r_avg/N << " " << r2/N << " " << -1+r_avg/N+r2/N-pow(r_avg/N,2.0) << " " << -1+2*avg_a*i*deltaT+avg_a2*i*i*deltaT*deltaT-avg_a*avg_a*i*i*deltaT*deltaT << " " << -1+2*avg_a*i*deltaT+avg_a2_1*i*i*deltaT*deltaT-avg_a*avg_a*i*i*deltaT*deltaT << " " << -1+2*avg_a*i*deltaT+avg_a2_2*i*i*deltaT*deltaT-avg_a*avg_a*i*i*deltaT*deltaT << endl;
        fyl.close();   
        cout << i << endl;
    }
 
    cout << r_avg << endl;    // pour afficher la valeur de r_avg à l'écran
}



