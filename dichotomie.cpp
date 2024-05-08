#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

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
	double delta = 100;
	double t_0 = 100.05;
	double t = 100;
	int newSign = -1;
	int sign = -1;
	
	while (delta >= 0.01)
	{
		sign = newSign;
		
		if (t - t_0 > 0)
	    {
	    	newSign = 1;
	    }
	    
	    else
	    {
	    	newSign = -1;
	    }
		
		if (newSign - sign != 0)
		{
			delta /= 10;
		}
		
		t += -newSign*delta;
		cout << t << " " << delta << endl;
	}  
	
    cout << t << " " << delta << endl;

   
}



