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
    double p = 1;
    double r1;
    double gam = 1.01;
    int correlations;
    string filename;
    
    
    
     
    cout << "system size?" << endl;
    cin >> N;
    cin.ignore();
    
    cout << "number of iterations?" << endl;
    cin >> Niter;
    cin.ignore();
    
    cout << "gamma" << endl;
    cin >> gam;
    cin.ignore();
    
    cout << "p? (1 for voter)" << endl;
    cin >> p;
    cin.ignore();
    
    cout << "1 for b=1, 2 for a and b indpdt, 3 for a=b, 4 for max. neg. corr., 5 for b=1/a" << endl;
    cin >> correlations;
    cin.ignore();
            
    cout << "enter filename.extension" << endl;
    getline(cin,filename); 
    
    cout << "hello world" << endl;
    
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
    
    ofstream fil (filename);
	fil.close();
	
	
	
	double epsPow1mGam = pow(eps,1-gam);
	int Seed = randdd.randInt();
	double tau(0);
	double tau_attempts = 0; 
	double sum_lambda_a = 0;
	double sum_lambda_b = 0;
	double avg_b = 0;
	double avg_a = 0;
	
	#pragma omp parallel for reduction(+:sum_lambda_a) reduction(+:sum_lambda_b) reduction(+:avg_a) reduction(+:avg_b) schedule(dynamic)
	for (int iter = 1; iter <= Niter; iter++)
	{
		MTRand mtrand1(iter+Seed);
	
		double r;
		
		
		vector <double> act(N);
	
		double sum_act = 0;
	
		vector <double> attrac(N);
		double sum_attrac = 0;
	
		double avg_a_sur_delta = 0;
		double avg_b_sur_delta = 0;
		double avg_ab_sur_delta = 0;
		double avg_a2_sur_delta = 0;
		double avg_b2_sur_delta = 0;
		


///////////////////////////////////////////////////////////////////   ordered list of random numbers between 0 and 1

		if (correlations == 4)
		{
			double list[N];
			
			for (int i = 0; i < N; i++)
			{
				r = mtrand1.randDblExc();
				list[i] = r;
			}
/*		
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
*/
			
			merge_sort(list,N);
			
			for (int i = 0; i < N; i++)
			{
				act[i] = pow(list[i]*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
				sum_act += act[i];					      
			}
		
			for (int i = 0; i < N; i++)
			{
				attrac[i] = act[N-1-i];      
				sum_attrac += attrac[i];					       
			}
		
			for (int i = 0; i < N; i++)
			{
				avg_b_sur_delta += attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a_sur_delta += act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_ab_sur_delta += attrac[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_b2_sur_delta += attrac[i]*attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a2_sur_delta += act[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);					       
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
				attrac[i] = 1;													
			}
		
			sum_attrac = N;
		
			for (int i = 0; i < N; i++)
			{
				avg_b_sur_delta += attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a_sur_delta += act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_ab_sur_delta += attrac[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_b2_sur_delta += attrac[i]*attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a2_sur_delta += act[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);					       
			}
		}
	
		if (correlations == 2)
		{			
			for (int i = 0; i < N; i++)
			{
				r = mtrand1.randDblExc();				
				act[i] = pow(r*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
				sum_act += act[i];					
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
			}
		
			for (int i = 0; i < N; i++)
			{
				avg_b_sur_delta += attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a_sur_delta += act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_ab_sur_delta += attrac[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_b2_sur_delta += attrac[i]*attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a2_sur_delta += act[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);					       
			}	
		}
	
		if (correlations == 3)
		{
			for (int i = 0; i < N; i++)
			{
				r = mtrand1.randDblExc();									
				act[i] = pow(r*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
				sum_act += act[i];					
		
				attrac[i] = act[i];            
				sum_attrac += attrac[i];					       
			}
		
			for (int i = 0; i < N; i++)
			{
				avg_b_sur_delta += attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a_sur_delta += act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_ab_sur_delta += attrac[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_b2_sur_delta += attrac[i]*attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a2_sur_delta += act[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);					       
			}
		}
		
		if (correlations == 5)
		{
			for (int i = 0; i < N; i++)
			{
				r = mtrand1.randDblExc();									
				act[i] = pow(r*(1-epsPow1mGam)+epsPow1mGam,1/(1-gam));
				sum_act += act[i];					
		
				attrac[i] = 1/act[i];            
				sum_attrac += attrac[i];					       
			}
		
			for (int i = 0; i < N; i++)
			{
				avg_b_sur_delta += attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a_sur_delta += act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_ab_sur_delta += attrac[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_b2_sur_delta += attrac[i]*attrac[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);
				avg_a2_sur_delta += act[i]*act[i]/(p*act[i]*sum_attrac+(1-p)*attrac[i]*sum_act);					       
			}
		}	
	
		//tau_attempts += (sum_act/N)*(sum_attrac/N)*pow(p*avg_b_sur_delta*avg_a2_sur_delta+(1-p*avg_ab_sur_delta)*avg_a_sur_delta,2)/avg_a2_sur_delta/(1+p*p*(avg_b2_sur_delta*avg_a2_sur_delta-pow(avg_ab_sur_delta,2)));
		//tau += (sum_attrac/N)*pow(p*avg_b_sur_delta*avg_a2_sur_delta+(1-p*avg_ab_sur_delta)*avg_a_sur_delta,2)/avg_a2_sur_delta/(1+p*p*(avg_b2_sur_delta*avg_a2_sur_delta-pow(avg_ab_sur_delta,2)));
		sum_lambda_a += avg_a2_sur_delta/(p*avg_b_sur_delta*avg_a2_sur_delta+(1-p*avg_ab_sur_delta)*avg_a_sur_delta);
		sum_lambda_b += (1-p*avg_ab_sur_delta)/(p*avg_b_sur_delta*avg_a2_sur_delta+(1-p*avg_ab_sur_delta)*avg_a_sur_delta);   //////////////  (1-p)*sum_lambda_b
		avg_a += sum_act/N;
		avg_b += sum_attrac/N;	
	}		
	
	cout << sum_lambda_a << " " << sum_lambda_b << " " << avg_a << " " << avg_b << endl;	
	
	//for (int i = 0; i < 50; i++)
	for (int i = 0; i < 100; i++)
	{   
/*
		double a = pow(eps,1-i/49.);		
		
		for (int j = 0; j < 50; j++)
		{
			double b = pow(eps,1-j/49.);
			ofstream fyl (filename,ios::app);
        	fyl << a << " " << b << " " << (p*b*sum_lambda_a+a*sum_lambda_b)/(p*a*avg_b+(1-p)*avg_a*b) << endl;
        	fyl.close();           
        }
*/
        double a_sur_b = 1/eps*pow(eps,2*(1-i/99.));
        
        ofstream fyol (filename,ios::app);
        fyol << a_sur_b << " " << (p*sum_lambda_a+a_sur_b*sum_lambda_b)/(p*a_sur_b*avg_b+(1-p)*avg_a) << endl;
        //fyol << endl;
        fyol.close();		
	}	
}


