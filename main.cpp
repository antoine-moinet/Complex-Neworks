#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "MersenneTwister.h"
#include <string>
//#include <boost/math/special_functions/erf.hpp>


using namespace std;

int main()
{   
    ofstream fil ("<k2>-<k>.dat");
    fil << 0 << " " << 0 << " " << 0 << endl;
    fil.close();

    ofstream files ("<r2>+3<r>2.dat");
    files << 0 << " " << 0 << endl;
    files.close();
    

    
    string filename;
    cout << "enter filename.extension" << endl;
    getline(cin,filename);

    
    MTRand mtrand1;

    long long int const N(10000000);
    int Niter = 1;
    double t(0.2);
    double deltaT = t/10;
    double ta(0);
    double tau;
    
    double alpha(0.9);
    double cmax = 1000000;
    double dmax(0);
    double r(0.5);
    //double r_avg(0);
    double betha(1.7);
    
    double sum(0);
    double Dmax(0);
    double c_0(1);
    
    double eps = 0.001;
    double gamma(1);
    //double r2(0);
    double A;
    double avg_act;
    double avg2_act;
    
    int const N1 = t/deltaT;

    double *act = new double[N];
    int *count = new int[N];
    double *Time = new double[N];
    //double *r2 = new double[N1];
    //double *r_avg = new double[N1];
    double r_avg;
    double r2;
   
    A = (betha)/(1-pow(cmax,-betha));
   
    //t = 1/(pow(avg2_act,0.5)+avg_act);
    
    vector<vector<double> > Tperc(20);  
    
    for (int u = 0; u <= 19; u++)
	{
		for (int v = 0; v <= 8; v++)
		{
		    Tperc[u].push_back(0);				    
		}
	}	
    
    Tperc[0][0] = 9.902790e-01; 
    Tperc[0][1] = 9.353840e-02; 
    Tperc[0][2] = 3.493670e-02; 
    Tperc[0][3] = 1.921520e-02; 
    Tperc[0][4] = 1.253570e-02; 
    Tperc[0][5] = 9.325020e-03; 
    Tperc[0][6] = 7.361060e-03; 
    Tperc[0][7] = 6.083520e-03; 
    Tperc[0][8] = 5.027710e-03; 
    Tperc[1][0] = 1.458700e+00; 
    Tperc[1][1] = 1.342940e-01; 
    Tperc[1][2] = 4.770920e-02; 
    Tperc[1][3] = 2.598030e-02; 
    Tperc[1][4] = 1.712410e-02; 
    Tperc[1][5] = 1.228870e-02; 
    Tperc[1][6] = 9.663050e-03; 
    Tperc[1][7] = 7.736550e-03; 
    Tperc[1][8] = 6.691880e-03; 
    Tperc[2][0] = 2.064860e+00; 
    Tperc[2][1] = 2.057840e-01; 
    Tperc[2][2] = 6.781840e-02; 
    Tperc[2][3] = 3.513840e-02; 
    Tperc[2][4] = 2.338460e-02; 
    Tperc[2][5] = 1.645600e-02; 
    Tperc[2][6] = 1.241160e-02; 
    Tperc[2][7] = 1.036010e-02; 
    Tperc[2][8] = 8.562070e-03; 
    Tperc[3][0] = 3.023160e+00; 
    Tperc[3][1] = 2.867570e-01; 
    Tperc[3][2] = 9.929290e-02; 
    Tperc[3][3] = 4.964640e-02; 
    Tperc[3][4] = 3.226340e-02; 
    Tperc[3][5] = 2.256660e-02; 
    Tperc[3][6] = 1.756920e-02; 
    Tperc[3][7] = 1.378930e-02; 
    Tperc[3][8] = 1.128330e-02; 
    Tperc[4][0] = 4.497710e+00; 
    Tperc[4][1] = 4.412570e-01; 
    Tperc[4][2] = 1.405540e-01; 
    Tperc[4][3] = 7.240620e-02; 
    Tperc[4][4] = 4.513310e-02; 
    Tperc[4][5] = 3.175050e-02; 
    Tperc[4][6] = 2.338460e-02; 
    Tperc[4][7] = 1.865000e-02; 
    Tperc[4][8] = 1.556240e-02; 
    Tperc[5][0] = 6.046320e+00; 
    Tperc[5][1] = 5.908930e-01; 
    Tperc[5][2] = 1.985860e-01; 
    Tperc[5][3] = 1.023020e-01; 
    Tperc[5][4] = 6.224980e-02; 
    Tperc[5][5] = 4.294260e-02; 
    Tperc[5][6] = 3.194400e-02; 
    Tperc[5][7] = 2.650250e-02; 
    Tperc[5][8] = 2.147130e-02; 
    Tperc[6][0] = 7.389210e+00; 
    Tperc[6][1] = 7.864780e-01; 
    Tperc[6][2] = 2.839180e-01; 
    Tperc[6][3] = 1.433790e-01; 
    Tperc[6][4] = 8.539220e-02; 
    Tperc[6][5] = 6.224980e-02; 
    Tperc[6][6] = 4.468630e-02; 
    Tperc[6][7] = 3.528600e-02; 
    Tperc[6][8] = 2.886410e-02; 
    Tperc[7][0] = 9.176230e+00; 
    Tperc[7][1] = 1.042750e+00; 
    Tperc[7][2] = 3.683220e-01; 
    Tperc[7][3] = 1.927450e-01; 
    Tperc[7][4] = 1.166110e-01; 
    Tperc[7][5] = 8.288080e-02; 
    Tperc[7][6] = 5.772820e-02; 
    Tperc[7][7] = 4.723690e-02; 
    Tperc[7][8] = 3.865220e-02; 
    Tperc[8][0] = 1.088450e+01; 
    Tperc[8][1] = 1.312960e+00; 
    Tperc[8][2] = 4.951390e-01; 
    Tperc[8][3] = 2.555520e-01; 
    Tperc[8][4] = 1.583790e-01; 
    Tperc[8][5] = 1.092220e-01; 
    Tperc[8][6] = 7.964680e-02; 
    Tperc[8][7] = 6.189240e-02; 
    Tperc[8][8] = 5.248020e-02; 
    Tperc[9][0] = 1.209260e+01; 
    Tperc[9][1] = 1.588690e+00; 
    Tperc[9][2] = 6.184350e-01; 
    Tperc[9][3] = 3.185870e-01; 
    Tperc[9][4] = 2.099210e-01; 
    Tperc[9][5] = 1.405540e-01; 
    Tperc[9][6] = 1.028920e-01; 
    Tperc[9][7] = 8.044330e-02; 
    Tperc[9][8] = 6.388800e-02; 
    Tperc[10][0] = 1.370500e+01; 
    Tperc[10][1] = 1.810900e+00; 
    Tperc[10][2] = 7.079010e-01; 
    Tperc[10][3] = 3.932390e-01; 
    Tperc[10][4] = 2.606880e-01; 
    Tperc[10][5] = 1.787450e-01; 
    Tperc[10][6] = 1.383630e-01; 
    Tperc[10][7] = 1.039210e-01; 
    Tperc[10][8] = 8.503490e-02; 
    Tperc[11][0] = 1.522620e+01; 
    Tperc[11][1] = 2.052350e+00; 
    Tperc[11][2] = 8.480790e-01; 
    Tperc[11][3] = 4.805770e-01; 
    Tperc[11][4] = 3.105170e-01; 
    Tperc[11][5] = 2.263630e-01; 
    Tperc[11][6] = 1.561550e-01; 
    Tperc[11][7] = 1.290540e-01; 
    Tperc[11][8] = 1.081410e-01; 
    Tperc[12][0] = 1.651890e+01; 
    Tperc[12][1] = 2.317000e+00; 
    Tperc[12][2] = 9.804740e-01; 
    Tperc[12][3] = 5.622140e-01; 
    Tperc[12][4] = 3.778950e-01; 
    Tperc[12][5] = 2.565440e-01; 
    Tperc[12][6] = 2.057840e-01; 
    Tperc[12][7] = 1.624950e-01; 
    Tperc[12][8] = 1.303440e-01; 
    Tperc[13][0] = 1.799080e+01; 
    Tperc[13][1] = 2.599930e+00; 
    Tperc[13][2] = 1.085090e+00; 
    Tperc[13][3] = 6.308660e-01; 
    Tperc[13][4] = 4.282800e-01; 
    Tperc[13][5] = 3.123100e-01; 
    Tperc[13][6] = 2.555520e-01; 
    Tperc[13][7] = 1.889470e-01; 
    Tperc[13][8] = 1.561550e-01; 
    Tperc[14][0] = 1.872130e+01; 
    Tperc[14][1] = 2.775810e+00; 
    Tperc[14][2] = 1.158500e+00; 
    Tperc[14][3] = 7.149800e-01; 
    Tperc[14][4] = 4.853830e-01; 
    Tperc[14][5] = 3.504460e-01; 
    Tperc[14][6] = 2.851110e-01; 
    Tperc[14][7] = 2.250630e-01; 
    Tperc[14][8] = 1.870770e-01; 
    Tperc[15][0] = 2.059350e+01; 
    Tperc[15][1] = 2.993230e+00; 
    Tperc[15][2] = 1.312960e+00; 
    Tperc[15][3] = 7.943430e-01; 
    Tperc[15][4] = 5.213770e-01; 
    Tperc[15][5] = 4.011430e-01; 
    Tperc[15][6] = 3.154330e-01; 
    Tperc[15][7] = 2.581070e-01; 
    Tperc[15][8] = 2.099210e-01; 
    Tperc[16][0] = 2.059350e+01; 
    Tperc[16][1] = 3.083930e+00; 
    Tperc[16][2] = 1.387910e+00; 
    Tperc[16][3] = 8.396830e-01; 
    Tperc[16][4] = 5.908930e-01; 
    Tperc[16][5] = 4.325630e-01; 
    Tperc[16][6] = 3.401390e-01; 
    Tperc[16][7] = 2.896250e-01; 
    Tperc[16][8] = 2.355540e-01; 
    Tperc[17][0] = 2.198660e+01; 
    Tperc[17][1] = 3.306390e+00; 
    Tperc[17][2] = 1.496610e+00; 
    Tperc[17][3] = 9.092560e-01; 
    Tperc[17][4] = 6.087970e-01; 
    Tperc[17][5] = 4.664440e-01; 
    Tperc[17][6] = 3.741530e-01; 
    Tperc[17][7] = 3.123100e-01; 
    Tperc[17][8] = 2.581070e-01; 
    Tperc[18][0] = 2.242850e+01; 
    Tperc[18][1] = 3.460500e+00; 
    Tperc[18][2] = 1.557380e+00; 
    Tperc[18][3] = 9.707660e-01; 
    Tperc[18][4] = 6.630460e-01; 
    Tperc[18][5] = 5.029780e-01; 
    Tperc[18][6] = 4.115690e-01; 
    Tperc[18][7] = 3.401390e-01; 
    Tperc[18][8] = 2.925210e-01; 
    Tperc[19][0] = 2.333920e+01; 
    Tperc[19][1] = 3.637020e+00; 
    Tperc[19][2] = 1.646280e+00; 
    Tperc[19][3] = 1.005960e+00; 
    Tperc[19][4] = 7.221300e-01; 
    Tperc[19][5] = 5.425470e-01; 
    Tperc[19][6] = 4.368890e-01; 
    Tperc[19][7] = 3.610650e-01; 
    Tperc[19][8] = 2.984010e-01; 


 

	
	
	ofstream filo (filename);
    //filo << -1 << " ";
    filo.close();   

    
    for (int u = 1; u <= 20; u++)
    {        
        //A = betha/(pow(eps,-betha)-1);
             

        for (int v = 1; v <= 9; v++)
    	{
            r_avg = 0;
            betha = 0.1*u;
            r2 = 0;
            alpha = 0.1*v;
            A = betha/(1-pow(cmax,-betha));
            double c2al_avg = A*(pow(cmax,1-betha)-1)/(1-betha); 
            //t = pow(pow(ta,0.5)*pow(tgamma(0.5),3)*tgamma(2.5)/2/c2al_avg,1/1.5);
            t = Tperc[u-1][v-1];
            cout << u << " " << v << endl;
            
            for (int k = 1; k <= Niter; k++)
            {  				   
			    for (int i = 0; i < N; i++)
				{
				    r = mtrand1.randDblExc();
				    //act[i] = pow(1-r,-1/betha);
				    act[i] = pow(1-r*betha/A,-1/betha);
				    //act[i] = pow(pow(eps,-betha)-r*betha/A,-1/betha);            
				    //act[i] = 1;       
				}			 
			        
		        for (int compteur = 0 ; compteur < N ; compteur++)
	    		{
		        	r = mtrand1.randDblExc();
		        	count[compteur] = 0;
		        	//tau = -log(1-r)/act[compteur];
		        	tau = (pow(1-r,-1/alpha)-1)/act[compteur];
		        	//tau = pow(boost::math::erf_inv(1-r),-2)/act[compteur];
		        	
		       	 	while (tau < t + ta)
		        	{
	            		if (tau > ta)
	            		{
	                		count[compteur] += 1;
	            		}                
	            		r = mtrand1.randDblExc();
	            		//tau += pow(boost::math::erf_inv(1-r),-2)/act[compteur];
	            		tau += (pow(1-r,-1/alpha)-1)/act[compteur];
	            		//tau += -log(1-r)/act[compteur];

		        	}

		        	r_avg += count[compteur]/(double)N;
		        	r2 += count[compteur]*count[compteur]/(double)N;
		    	}		    
            }
            
            //double gam = pow(tgamma(1-alpha),2.)*tgamma(1+alpha)*tgamma(1+2*alpha);
            
            ofstream fyl (filename,ios::app);
        	//fyl << betha << " " << alpha << " " << t_perc/pow(3*gam/(2*c2_avg*tgamma(1+alpha)-c_avg*c_avg),1/(2*alpha)) << endl;
        	//fyl << betha << " " << alpha << " " << t_perc*alpha*2*c_avg << endl;
        	//fyl << betha << " " << alpha << " " << 1-r_avg-r2+r_avg*r_avg << endl;
            fyl << betha << " " << alpha << " " << r_avg << endl;
        	fyl.close();           
        }

        ofstream fyol (filename,ios::app);
        fyol << endl;
        fyol.close();
    }

/*    
    for (int v = 1; v <= t/deltaT; v++)
	{
		r_avg[v-1] = 0;
		r2[v-1] = 0;
	}
    
    for (int k = 1; k <= Niter; k++)  
    {
    	for (int i = 0; i < N; i++)
		{
		    r = mtrand1.randDblExc();
		    act[i] = c_0*pow(1-r,-1/betha);
		    //act[i] = c_0;
		    //act[i] = pow(pow(eps,-betha)-r*betha/A,-1/betha);
		    //act[i] = pow(1-r*betha/A,-1/betha);
		    r = mtrand1.randDblExc();
		    Time[i] = (pow(1-r,-1/alpha)-1)/(gamma*act[i]);
		    count[i] = 0;		    
		}
		
		//cout << k << endl;
		
		
		
		//int v = 0;
		for (int v = 1; v <= t/deltaT; v++)
		//while (r_avg+r2 < 1)
		{
		   
		    //v++;
		    
		    cout << v << endl;
		    
		    for (int compteur = 0 ; compteur < N ; compteur++)
		    {
		        r = mtrand1.randDblExc();
		        //cout << act[compteur] << endl;
		        
		        //tau = -log(1-r)/act[compteur];
		        //tau = (pow(1-r,-1/alpha)-1)/(gamma*act[compteur]);
		        //tau = pow(boost::math::erf_inv(1-r),-2)/act[compteur];
		        while (Time[compteur] < ta + v*deltaT)
		        {
		            if (Time[compteur] > ta)
		            {
		                count[compteur] += 1;                                     
		            }

		            r = mtrand1.randDblExc();   
		            Time[compteur] += (pow(1-r,-1/alpha)-1)/(act[compteur]*gamma);
		            //tau += pow(boost::math::erf_inv(1-r),-2)/act[compteur];
		            //tau += -log(1-r)/act[compteur];
		        }

		        r_avg[v-1] += count[compteur]/(double)N;
		        r2[v-1] += count[compteur]*count[compteur]/(double)N;
		    }		    
		} 
	} 
        
*/    

        //cout << r_avg << " " << r2-3*r_avg+4*r_avg*r_avg << " " << r2-1+r_avg << endl;
/*
        //cout << "liste des nombres de tirs générée" << endl;

        for (int compteur = 0 ; compteur < N ; compteur++)
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
*/     
/*
        for (int i = 0 ; i < d ; i++)
        {
            sum += pr[i];
        }        
        
        cout << v << endl;
        
        for (int i = 0 ; i < d ; i++)
        {
            r_avg += i*pr[i];
        }

        double r2(0);

        for (int i = 0 ; i < d ; i++)
        {
            r2 += i*i*pr[i];
        }
*/        
/*
        ofstream fil ("<k2>-<k>.dat",ios::app);
        //fil << v*deltaT << " " << 3*r_avg*r_avg+r2 << endl;
        fil << v*deltaT << " " << r_avg+r2 << " " << r2 << endl;
        fil.close();
*/
/*
    ofstream  fyl ("1-<r2>-<r>.dat");
	for (int v = 1; v <= t/deltaT; v++)
	{
		//fyl << v*deltaT << " " << r2[v-1]/(double)Niter << " " << 1-(r_avg[v-1]+r2[v-1]+r_avg[v-1]*r_avg[v-1]/(double)Niter)/(double)Niter << " " << r_avg[v-1]/(double)Niter << endl;
		fyl << v*deltaT << " " << r2[v-1]-r_avg[v-1]*r_avg[v-1] << " " << r_avg[v-1] << endl;
	}
    fyl.close();
        
 */       
       

        
            
/*        
        
        double *G = new double[N];

        for (int k = 0; k < N; k++)
        {
            G[k] = 0;
        }
            
        for (int k = 0; k < N; k++)
        {
            while (count[k] != 0)
            {
                int n = mtrand1.randInt(N-1);
                while (n == k)
                {
                    n = mtrand1.randInt(N-1);
                }

                G[k] += 1;
                G[n] += 1;
                count[k] -= 1;
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


        int Nbin = 50;
        double a = pow(10,5/(double)Nbin);        
        double *g = new double[Nbin];


        for (int k = 0; k < Nbin; k++)
        {
            g[k] = 0;
            double count = 1/(double)N;
            for (int i = (int)(10*pow(a,k)); i < (int)(10*pow(a,k+1)); i++)
            {
                count += 1;
                if (i < D)
                {
                    g[k] += p[i];
                }
                
            }
            g[k] = g[k]/count;
            //cout << g[k] << " " << k << endl;
        }

*/        
/*     
        for (int i = 0 ; i < D ; i++)
        {
            sum += p[i];
        }

        //cout << sum << endl;
        

        double k(0);
        double k2(0);

        for (int i = 0 ; i < D ; i++)
        {
            k2 += p[i]*i*(i-1);
            k += p[i]*i;
        }

        ofstream files ("<k2>-<k>.dat",ios::app);
        files << v*deltaT << " " << k << " " << k2 << endl;
        files.close();   




        

        

*/
/*
        ofstream out ("Pr.dat");

        for (int compteur = 0 ; compteur < d; compteur++)
        {
            out << compteur << " " << pr[compteur] << " " << erf((compteur+1)/sqrt(t))-erf(compteur/sqrt(t)) << endl;
        }
        out.close();


        ofstream op ("BinDegDistr.dat");        

        for (int compteur = 0 ; compteur < Nbin; compteur++)
        {
            op << 10*pow(a,compteur)-r_avg << " " << g[compteur] << " " << endl;
        }
        op.close();


        ofstream pop("degDistr.dat");

        for (int compteur = 0 ; compteur < D; compteur++)
        {
            pop << compteur << " " << p[compteur] << " " << endl;
        }
        pop.close();


        delete [] p;
        delete [] G;
        //delete [] g;
        //delete [] pr;
*/
    
    
    delete [] count;    
    delete [] act; 
    delete [] Time;
    //cout << r_avg << " " << sin(3.14*alpha)/(3.14*alpha)*pow(avg_act*t,alpha) << endl;  
}



