#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "MersenneTwister.h"
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
   

    ofstream fil ("perctime.dat");
    fil.close();

    ofstream  fyl ("<r2>+3<r>2-<r>.dat");
    fyl << 0 << " " << 0 << " " << 1 << endl;
    fyl.close();
    
    MTRand mtrand1;

    int const N(1000000);
    double t(0.4);
    //double tau;
    double ta(0);
    
    double alpha(0.5);
    double dmax(0);
    double r(0.5);
    double r_avg(0);
    double deltaT(0.1);
    double sum(0);
    double Dmax(0);
    double eps(0.1);
    double betha(0.01);
    double A;
    double gamma;
    int x;

    double *act = new double[N];
    int *count = new int[N];
    double *tau = new double[N];
    int *countTot = new int[N];

    double r2(0);

    ofstream filo ("perctime.dat",ios::app);
    filo << 0 << " ";
    filo.close();

    gamma = 1.001;

    for (int i = 0; i < N; i++)
    {
        r = mtrand1.randDblExc();
        //act[i] = pow(1-r,-1/betha);
        act[i] = pow((1-pow(eps,1-gamma))*r+pow(eps,1-gamma),1/(1-gamma));
        //act[i] = (1-eps)*r+eps;
        //act[i] = (1/eps-1)*r+1;
        //act[i] = 1;  
        countTot[i] = 0;    
    }
    
    

    for (int v = 1; v <= 50; v++)
    {
        r_avg = 0;
        r2 = 0;
        sum = 0;
        dmax = 0;
        Dmax = 0;
        x = 1;
        alpha = 0.02*v+0.001;

        for (int i = 0; i < N; i++)
        {
            r = mtrand1.randDblExc();
            tau[i] = (pow(1-r,-1/alpha)-1)/act[i]; 
            countTot[i] = 0;   
        }
        
        while (1-r_avg-r2 > 0)
        {
            for (int compteur = 0 ; compteur < N ; compteur++)
            {
                r = mtrand1.randDblExc();
                count[compteur] = 0;
                //tau = -log(1-r)/act[compteur];
                //tau[compteur] += (pow(1-r,-1/alpha)-1)/act[compteur];
                while (tau[compteur] < x*deltaT + ta)
                {
                    if (tau[compteur] > ta)
                    {
                        count[compteur] += 1;
                        countTot[compteur] += 1;
                    }                
                    r = mtrand1.randDblExc();
                    tau[compteur] += (pow(1-r,-1/alpha)-1)/act[compteur];
                    //tau += -log(1-r)/act[compteur];

                }

                r_avg += count[compteur]/(double)N;
            }

            r2 = 0;

            for (int i = 0 ; i < N ; i++)
            {
                r2 += pow(countTot[i]-r_avg,2)/(double)N;
            }

            x++;
            //cout << 1-r_avg-r2 << " " << countTot[0] << " " << r_avg << " " << r2 << endl;
        }        

        

        ofstream fyl ("perctime.dat",ios::app);
        fyl << x*deltaT << " ";
        fyl.close();
    }

    ofstream fyol ("perctime.dat",ios::app);
    fyol << endl;
    fyol.close();

    for (int u = 1; u <= 100; u++)
    {
        gamma = 1 + u*0.02 + 0.001;

        ofstream fil ("perctime.dat",ios::app);
        fil << gamma << " ";
        fil.close();

        for (int i = 0; i < N; i++)
        {
            r = mtrand1.randDblExc();
            //act[i] = pow(1-r,-1/betha);
            act[i] = pow((1-pow(eps,1-gamma))*r+pow(eps,1-gamma),1/(1-gamma));
            //act[i] = pow(1-r*(0.1*u+0.01)/A,-1/(0.1*u+0.01));
            //act[i] = (1-eps)*r +  eps;;       
        }        

        for (int v = 1; v <= 50; v++)
        {
            r_avg = 0;
            r2 = 0;
            sum = 0;
            dmax = 0;
            Dmax = 0;
            x = 1; 
            alpha = 0.02*v+0.001;

            for (int i = 0; i < N; i++)
            {
                r = mtrand1.randDblExc();
                tau[i] = (pow(1-r,-1/alpha)-1)/act[i];    
                countTot[i] = 0;    
            }
            
            while (1-r_avg-r2 > 0)
            {
                for (int compteur = 0 ; compteur < N ; compteur++)
                {
                    r = mtrand1.randDblExc();
                    count[compteur] = 0;
                    //tau = -log(1-r)/act[compteur];
                    //tau[compteur] += (pow(1-r,-1/alpha)-1)/act[compteur];
                    while (tau[compteur] < x*deltaT + ta)
                    {
                        if (tau[compteur] > ta)
                        {
                            count[compteur] += 1;
                            countTot[compteur] += 1;
                        }                
                        r = mtrand1.randDblExc();
                        tau[compteur] += (pow(1-r,-1/alpha)-1)/act[compteur];
                        //tau += -log(1-r)/act[compteur];
                    }

                    r_avg += count[compteur]/(double)N;
                }    

                r2 = 0;            

                for (int i = 0 ; i < N ; i++)
                {
                    r2 += pow(countTot[i]-r_avg,2)/(double)N;
                }

                x++;
                //cout << x << endl;
            }


            

            ofstream fyl ("perctime.dat",ios::app);
            fyl << x*deltaT << " ";
            fyl.close();
        }

        ofstream fyol ("perctime.dat",ios::app);
        fyol << endl;
        fyol.close();
    }

    
    
  
       
        

        //cout << r_avg << endl;

        //cout << "liste des nombres de tirs générée" << endl;

/*        for (int compteur = 0 ; compteur < N ; compteur++)
        {
            if (dmax < count[compteur])
            {
                dmax = count[compteur];
            }
        }
        
        
        int const d(dmax+1);
        double *pr = new double[d];

        for (int i = 0 ; i < d ; i++)
        {
            pr[i] = 0;
        }

        for (int i = 0 ; i < N ; i++)
        {
            int j = 0;
            while (count[i] != j)
            {
                j++;
            }
            pr[j] += 1/(double)N;
        }

        for (int i = 0 ; i < d ; i++)
        {
            sum += pr[i];
        }  

        r_avg = 0;
        for (int i = 0 ; i < d ; i++)
        {
            r_avg += i*pr[i];
        }

        r2 = 0;

        for (int i = 0 ; i < d ; i++)
        {
            r2 += i*i*pr[i];
        }
       
        ofstream fil ("<r2>+3<r>2.dat",ios::app);
        //fil << v*deltaT << " " << 3*r_avg*r_avg+r2 << endl;
        fil << v*deltaT << " " << r_avg << " " << r2+r_avg*r_avg << endl;
        fil.close();

        ofstream fyl ("<r2>+3<r>2-<r>.dat",ios::app);
        fyl << v*deltaT << " " << 3*r_avg*r_avg+r2-3*r_avg << " " << r_avg*r_avg-r2-r_avg+1 <<  endl;
        fyl.close();
*/

/*        
        int *G = new int[N];

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

        

        for (int k = 0; k < N; k++)
        {
            if (G[k] > Dmax)
            {
                Dmax = G[k];
            }
        }

        int const D(Dmax+1);
        double *p = new double[D];


        for (int i = 0 ; i < D ; i++)
        {
            p[i] = 0;
        }

        for (int i = 0; i < N; i++)
        {
            int j = 0;
            while (G[i] != j)
            {
                j++;
            }
            p[j] += 1/(double)N;
        }

*/       

/*        double *pth = new double[D];
        sum = 0;

        for (int k = 0; k < D; k++)
        {
            pth[k] = 0;
        }

        for (int k = 0; k < D; k++)
        {
            for (int i = 0; i <= k; i++)
            {
                pth[k] += pr[i]*exp(-r_avg+(k-i)*log(r_avg)-logsum(k-i));                
            }
        }     

        double kmoy(0);
        double k2(0);

        for (int i = 0 ; i < D ; i++)
        {
            k2 += p[i]*i*(i-1);
            kmoy += p[i]*i;
        }

        double *g = new double[D];
        int Nr = 0;
        for (int k = 0; k < D; k++)
        {
            g[k] = 0;
        }

        for (int i = 0; i < N; i++)
        {
            if (count[i] == int(r_avg))
            {
                g[G[i]] += 1;
                Nr += 1;
            }
        }

        sum = 0;

        for (int k = 0; k < D; k++)
        {
            g[k] /= (double)Nr;
            sum += k*g[k];
        }



        ofstream out ("Pr.dat");

        for (int compteur = 0 ; compteur < d; compteur++)
        {
            out << compteur << " " << pr[compteur] << " " << pow(t,-0.5)*exp(-atan(1)*compteur*compteur/t) << endl;
        }
        out.close();
*/
/*
        ofstream op ("AgeDegDistr.dat");

        for (int compteur = 0 ; compteur < D; compteur++)
        {
            op << compteur << " " << p[compteur] << endl;
        }
        op.close();





        double r_avg2(0);

        for (int compteur = 0 ; compteur < N ; compteur++)
        {
            r = mtrand1.randDblExc();
            count[compteur] = 0;
            //tau = -log(1-r)/act[compteur];
            tau = (pow(1-r,-1/alpha)-1)/act[compteur];
            while (tau < v*deltaT)
            {            
                count[compteur] += 1;                
                r = mtrand1.randDblExc();
                tau += (pow(1-r,-1/alpha)-1)/act[compteur];
                //tau += -log(1-r)/act[compteur];
            }

            r_avg2 += count[compteur]/(double)N;
        }

        double Dravg = r_avg2-r_avg;

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
        
        int const DD(Dmax+1);
        double *pp = new double[DD];

        for (int i = 0 ; i < DD ; i++)
        {
            pp[i] = 0;
        }

        for (int i = 0; i < N; i++)
        {
            int j = 0;
            while (G[i] != j)
            {
                j++;
            }
            pp[j] += 1/(double)N;
        }

        cout <<  Dravg << endl;
        cout << r_avg << endl;

        double gam = 1.165;    //gamma(1-alpha)*gamma(1+alpha)

        ofstream opa ("degDistr.dat");

        for (int k = 0 ; k < DD; k++)
        {
            //opa << k << " " << pp[k] << " " << pp[k+(int)Dravg] - pow(ta,alpha)/gam*(pp[k+(int)Dravg-1]-pp[k+(int)Dravg]) << endl;
            opa << k << " " << pp[k] << " " << exp(-r_avg+k*log(r_avg)-logsum(k))  << endl;
        }
        opa.close();


        
        //delete[] pth;
        //delete[] g;
        delete [] p;
        delete [] G;
        delete [] pp;
        //delete [] pr;
*/        
    

    
    delete [] count;    
    delete [] act; 
    delete [] tau; 
    delete [] countTot; 
}



