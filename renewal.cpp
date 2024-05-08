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
    files << 0 << " " << 0 << endl;
    files.close();

    
    MTRand mtrand1;          // ça c'est pour pouvoir utiliser twister il faut le déclarer

    int const N(1000000);  // tu déclare toutes tes variables en leur assignant ou non une valeur les () ou = c'est pareil
    int N1 = 100;
    
    double deltaT = 0.01;
    double r_avg, tau, r;
    double ta = 0;
    double alpha = 0.9;
    double act = 5;
    tau = 0;
    int renewals = 0;
    
    while (tau < 1000)
    {
    	r = mtrand1.randDblExc();
        //tau = pow(boost::math::erf_inv(1-r),-2);
        tau += (pow(1-r,-1/alpha)-1)/act;
        renewals++;
        
        ofstream fyl ("r_avg.dat",ios::app);
        fyl << tau << " " << renewals << endl;
        fyl.close(); 
    }
    
    ofstream files1 ("r_avg1.dat");
    files1 << 0 << " " << 0 << endl;
    files1.close();
    
    renewals = 0;
    tau = 0;
    
    while (tau < 1000)
    {
    	r = mtrand1.randDblExc();
        //tau = pow(boost::math::erf_inv(1-r),-2);
        tau += -log(1-r);
        renewals++;
        
        ofstream fyl1 ("r_avg1.dat",ios::app);
        fyl1 << tau << " " << renewals << endl;
        fyl1.close(); 
    }
    
    ofstream files2 ("r_avg2.dat");
    files2 << 0 << " " << 0 << endl;
    files2.close();
    
    renewals = 0;
    tau = 0;
    
    while (tau < 1000)
    {
    	r = mtrand1.randDblExc();
        //tau = pow(boost::math::erf_inv(1-r),-2);
        tau += 1;
        renewals++;
        
        ofstream fyl2 ("r_avg2.dat",ios::app);
        fyl2 << tau << " " << renewals << endl;
        fyl2.close(); 
    }
    
}



