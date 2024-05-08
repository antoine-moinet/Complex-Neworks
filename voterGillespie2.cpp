#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include "MersenneTwister.h"
using namespace std;


void merge(double a[], int left_low, int left_high, int right_low, int right_high) 
{ 
    int length = right_high-left_low+1;
    double* temp = new double [length];
    int left = left_low;
    int right = right_low;
    for (int i = 0; i < length; ++i) { 
        if (left > left_high)
            temp[i] = a[right++];
        else if (right > right_high)
            temp[i] = a[left++];
        else if (a[left] <= a[right])
            temp[i] = a[left++];
        else
            temp[i] = a[right++]; 
    }
    
    for (int i=0; i< length; ++i) 
        a[left_low++] = temp[i];
        
    delete [] temp;
}


void Merge_sort(double a[], int low, int high) {
    if (low >= high)                  //Base case: 1 value to sort->sorted
      return;                         //(0 possible only on initial call)
    else {
      int mid = (low + high)/2;       //Approximate midpoint*
      Merge_sort(a, low, mid);        //Sort low to mid part of array
      Merge_sort(a, mid+1, high);     //Sort mid+1 to high part of array
      merge(a, low, mid, mid+1,high); //Merge sorted subparts of array
    }
}

void merge_sort(double a[], int length) {
    Merge_sort(a, 0, length-1);
}






int main()
{  
    MTRand randdd;
	
    

    int N(100);   //number of vertices
    int Niter = 10000;
    double par = 1;
    double r1;
    int correlations;
    string filename;
    
    
    
     
    cout << "system size?" << endl;
    cin >> N;
    cin.ignore();
    
    cout << "number of iterations?" << endl;
    cin >> Niter;
    cin.ignore();
    
    cout << "1 for b=1, 2 for a and b indpdt, 3 for b=a^2, 4 for max. neg. corr., 5 for b=1/a" << endl;
    cin >> correlations;
    cin.ignore();
            
    cout << "enter filename.extension" << endl;
    getline(cin,filename); 
    
/*    
    double list[N];
    double r;
					
	for (int i = 0; i < N; i++)
	{
		r = randdd.randDblExc();
		list[i] = r;
	}

	
	merge_sort(list,N);
*/       
    
    double eps = 0.001;
    
    
    vector<int> permutation(N);
    int local_count = 1;
    int i = 0;
/*    
    while (i < N-local_count)
    {
    	r1 = randdd.randDblExc();
    	if (r1 < 0.5)
    	{
    		permutation[i] = N - local_count;
    		permutation[N-local_count] = i;
    		local_count++;   		
    	}
    	
    	else 
    	{
    		permutation[i] = i+1;
    		permutation[i+1] = i;
    		i++;
    	}
    	i++;
    }
*/

	for (int i = 0; i < N; i++)
    {
    	permutation[i] = i;    	
    }
    
	vector<int> vu(N,0);
	int vus = 0;	    
    while (vus != N)
    {
    	while (vu[i] == 1)
    	{
    		i++;
    	}
    	
    	local_count = randdd.randInt(N-1-i-1)+i+1;
    	
    	while (vu[local_count] == 1)
    	{
    		local_count = randdd.randInt(N-1-i-1)+i+1;
    	}
    	
    	permutation[i] = local_count;
    	permutation[local_count] = i;
        vu[local_count] = 1;
        vu[i] = 1;
        vus += 2;
        
    }
    
    local_count = 0;
		    
    for (int i = 0; i < N; i++)
    {
    	
    	if (permutation[permutation[i]] != i)
    	{
    		local_count++;
    	}
    }
    
    cout << local_count << " problems" << endl;
    
    ofstream fil (filename);
	fil.close();
	
	for (int f = 0; f <= 5; f++)
	{
		par = 0.2*f;
		
		for (int i = 0; i <= 40; i++)
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
						   
				for (int i = 0; i < N; i++)
				{              
				    r = mtrand1.randDblExc();
					if (r < 0.5)
					{
						opinion[i] = 1;
						sum++;
					}
				    
				}                   

				//while(t < N*N*50)
				while (sum != 0 && sum != N)
				{
				
		////////////////////////////////////////////////////// choose shooter     
			   
				    r = mtrand1.randDblExc();
				    int j = r*(N-1);
				    
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
				    int Nguy = r*(N-1);
				   
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
				    t_steps += -log(1-r)/N; 
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
			files << par << " " << gam << " " << t_steps/Niter/N/log(2) << endl;
		
			files.close();
			cout << par << " " << Exit/(double)Niter << " " << t_steps/Niter/N/log(2) << " Niter = " << Niter << " N = " << N << endl;
		
		}

        ofstream fyol (filename,ios::app);
        fyol << endl;
        fyol.close();		
	}
	
}


