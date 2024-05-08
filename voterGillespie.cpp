#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include "MersenneTwister.h"
using namespace std;



int main()
{  
    MTRand randdd;
	
    

    int N(100);   //number of vertices
    int Niter = 10000;
    double par = 1;
    double r1;
    int correlations;
    string filename;
    double eps = 0.001;
     
    cout << "system size?" << endl;
    cin >> N;
    cin.ignore();
    
    cout << "number of iterations?" << endl;
    cin >> Niter;
    cin.ignore();
    
    cout << "1 for b=1, 2 for a and b indpdt, 3 for a=b, 4 for max. neg. corr." << endl;
    cin >> correlations;
    cin.ignore();
        
    cout << "p? (p=1 for voter)" << endl;
    cin >> par;
    cin.ignore();
    
    cout << "epsilon?" << endl;
    cin >> eps;
    cin.ignore();
            
    cout << "enter filename.extension" << endl;
    getline(cin,filename);
    
    
   
    
    
    
    
    
    
    
    ofstream fil (filename);
	fil.close();
	
/*
	vector<int> magnetization(N*N*50,0);
	vector<int> m_1(N*N*50,0);
	vector<int> m_2(N*N*50,0);
	vector<int> m_3(N*N*50,0);
	vector<int> activeBonds(N*N*50,0);
	vector<int> survive(N*N*50,0);
	
   
   
    vector <double> act(N);
    vector <double> cum_act(N);
    double sum_act = 0;
    vector <double> cum_attrac(N);
    vector <double> attrac(N);
    double sum_attrac = 0;
    double avg_oneOverAct = 0;
   
    
    for (int i = 0; i < 50; i++)
    {
        act[i] = 1;      
    }
    for (int i = 50; i < 85; i++)
    {
        act[i] = 2;      
    }  
    for (int i = 85; i < 100; i++)
    {
        act[i] = 3;      
    }
    
    for (int i = 0; i < N; i++)
    {
        r1 = randdd.randDblExc();
        
        //act[i] = pow(r1*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
        
        attrac[i] = 4;
        sum_act += act[i];
        cum_act[i] = sum_act;
        sum_attrac += attrac[i];
        cum_attrac[i] = sum_attrac;  
        avg_oneOverAct += 1/act[i];      
    }    
    
    double omega = 0;
    vector<int> opinion(N,0);
    for (int i = 0; i < N; i++)
    {              
        if (act[i] < 2*eps)
        {
        	opinion[i] = 1;
       		omega += 1/act[i];                   
        }
         
    }                   
    cout << omega/avg_oneOverAct << endl;
*/

	for (int i = 0; i <= 35; i++)
	{
		double gam = 0.01 + i*0.1;
		double epsPow1mGam = pow(eps,1-gam);
		int Seed = randdd.randInt();
		double t(0); 
		double t_steps = 0; 
		int Exit = 0; 
		
		#pragma omp parallel for reduction(+:Exit) reduction(+:t) reduction(+:t_steps) schedule(dynamic)
		for (int iter = 1; iter <= Niter; iter++)
		{
		    MTRand mtrand1(iter+Seed);
			
		    double r;
		    int sum = 0;
		    //int sum_1 = 0;
		    //int sum_2 = 0;
		    //int sum_3 = 0;
		    
		    
		    
		    vector <double> act(N);
			vector <double> cum_act(N);
			double sum_act = 0;
			vector <double> cum_attrac(N);
			vector <double> attrac(N);
			double sum_attrac = 0;
		
		
	///////////////////////////////////////////////////////////////////   ordered list of random numbers between 0 and 1

			if (correlations == 4)
			{
				vector<double> list(N);
			    
				for (int i = 0; i < N; i++)
				{
					r = mtrand1.randDblExc();
					list[i] = r;
				}
		
				int a;
				double b;
				for (int k = 1; k < list.size() ; k++)
				{   
					a = k;
					while (a > 0 && list[a] < list[a-1])
					{            
						b = list[a];  
						//j = qui[u];              
						list[a] = list[a-1]; 
						//qui[u] = qui[u-1];               
						list[a-1] = b;
						//qui[u-1] = j;
						a -= 1;                           
					} 
					//cout << k << endl;           
				}
				
				for (int i = 0; i < N; i++)
				{
					act[i] = pow(list[i]*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
					sum_act += act[i];
					cum_act[i] = sum_act;       
				}
				
				for (int i = 0; i < N; i++)
				{
					attrac[i] = act[N-1-i];      
					sum_attrac += attrac[i];
					cum_attrac[i] = sum_attrac;        
				}
			}		
		
	//////////////////////////////////////////////////////////////////// 
		   
			if (correlations == 1)
			{
				for (int i = 0; i < N; i++)
				{
					r = mtrand1.randDblExc();
					act[i] = pow(r*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
					sum_act += act[i];
					cum_act[i] = sum_act;								
				}
			}
			
			if (correlations == 2)
			{			
				for (int i = 0; i < N; i++)
				{
					r = mtrand1.randDblExc();				
					act[i] = pow(r*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
					sum_act += act[i];
					cum_act[i] = sum_act;	
					
					attrac[i] = act[i];							     
				}
				
				for (int i = 1; i < N; i++)///////////////////////////////////// mélange de Fisher-Yates
				{
					int j =	mtrand1.randInt(N-i);
					double a_of_j = attrac[j];
					attrac[j] = attrac[N-i];
					attrac[N-i] = a_of_j;					     
				}			
				
				for (int i = 0; i < N; i++)
				{
					sum_attrac += attrac[i];
					cum_attrac[i] = sum_attrac;
				}	
			}
			
			if (correlations == 3)
			{
				for (int i = 0; i < N; i++)
				{
					r = mtrand1.randDblExc();									
					act[i] = pow(r*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
					sum_act += act[i];
					cum_act[i] = sum_act;
				
					attrac[i] = act[i];            
					sum_attrac += attrac[i];
					cum_attrac[i] = sum_attrac;        
				}
			}			
		    
		    vector<int> opinion(N,0);
		    double a_sur_b_max = 0;
		    double b_sur_a_moy = 0;
		    int myMan = 0;
		    		   
		    for (int i = 0; i < N; i++)
		    {              
		        //r = mtrand1.randDblExc();
				//if (r < 0.5)
				//{
				//	opinion[i] = 1;
				//	sum++;
				//}
				if (act[i]/attrac[i] > a_sur_b_max)
				{
					a_sur_b_max = act[i]/attrac[i];
					myMan = i;
				}
		        b_sur_a_moy += attrac[i]/act[i];
		    } 
		    
		    double omega = attrac[i]/act[i]/b_sur_a_moy;
		    //cout << omega << endl;
		    opinion[myMan] = 1;
		    sum = 1;
		                      

		    //while(t < N*N*50)
		    while (sum != 0 && sum != N)
		    {
		    
	////////////////////////////////////////////////////// choose shooter     
		   
		        r = mtrand1.randDblExc();
		        int j = r*N;
		        
		        r = r*sum_act;

				if (r > cum_act[j])
				{
					while (r > cum_act[j])
					{
						j++;								
					}
				}

				else
				{							
					while (r < cum_act[j-1] && j > 0)
					{
						j -= 1;
					}							
				}
		
	///////////////////////////////////////////////////////
		        
	/////////////////////////////////////////////////////// choose receiver
		        
		        r = mtrand1.randDblExc();
		        int Nguy = r*N;
		       
		        if (correlations != 1)
		        {
		        	r = r*sum_attrac;		       

					if (r > cum_attrac[Nguy])
					{
						while (r > cum_attrac[Nguy])
						{
							Nguy++;								
						}
					}

					else
					{							
						while (r < cum_attrac[Nguy-1] && Nguy > 0)
						{
							Nguy -= 1;
						}							
					}
		        }
		        
			
	////////////////////////////////////////////////////// 
			
		        r = mtrand1.randDblExc();
		        
		        if (r < par)
		        {
		            sum += opinion[Nguy]-opinion[j];
		            //sum_1 += (opinion[Nguy]-opinion[j])*(act[j]-2)*(act[j]-3)/2;
		            //sum_2 += -(opinion[Nguy]-opinion[j])*(act[j]-1)*(act[j]-3);
		            //sum_3 += (opinion[Nguy]-opinion[j])*(act[j]-1)*(act[j]-2)/2;
		            opinion[j] = opinion[Nguy];      // par = 1 voter model
		        }
		        else
		        {
		            sum += opinion[j]-opinion[Nguy];
		            //sum_1 += -(opinion[Nguy]-opinion[j])*(act[Nguy]-2)*(act[Nguy]-3)/2;
		            //sum_2 += (opinion[Nguy]-opinion[j])*(act[Nguy]-1)*(act[Nguy]-3);
		            //sum_3 += -(opinion[Nguy]-opinion[j])*(act[Nguy]-1)*(act[Nguy]-2)/2;
		            opinion[Nguy] = opinion[j];      // par = 0 invasion process
		        }  

		        r = mtrand1.randDblExc();               
		        t += -log(1-r)/sum_act; 
		        t_steps += -log(1-r)/N/((1-omega)*log(1/(1-omega))+omega*log(1/omega)); 
		        //t++; 
		        //t += 1./(N1*mu1);  
		       
		        //magnetization[t] += sum;
		        //m_1[t] += sum_1;
		        //m_2[t] += sum_2;
				//m_3[t] += sum_3;

					        
		        //activeBonds[t] += sum*(N-sum); 
		        
		        //if (sum != 0 && sum != N)
		        //{
		        //	survive[t]++; 
		        //}                                     
		        
		          
		    }

		    if (sum == N)
		    {
		        Exit++;
		    }
		}
		
		
		ofstream files (filename,ios::app);
		files << gam << " " << t/Niter << " " << t_steps/Niter/N << " " << Exit/(double)Niter << endl;
		
		files.close();
		cout << gam << " " << t/Niter << " " << t_steps/Niter/N << " " << Exit/(double)Niter << endl;
		
	}

    
    
    //cout << Exit/(double)Niter << " " << t/Niter << " Niter = " << Niter << " N = " << N << endl;
/*    
    int lastSurvive = 1;
    while (survive[lastSurvive] != 0)
    {
    	lastSurvive++;
    }
    
   
    ofstream files ("VoterOneRun.dat");
    for (int k = 1; k < lastSurvive; k++)
    {
    	files << k << " " << 2*activeBonds[k]/(double)(survive[k]*N*N) << " " << magnetization[k]/(double)(Niter*N) << " " << m_1[k]/(double)(Niter*50) << " " << m_2[k]/(double)(Niter*35) << " " << m_3[k]/(double)(Niter*15) << endl;
    }
    
    cout << survive[1] << endl;
    cout << survive[100] << endl;
    
    files.close();
*/  
      
}


