#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "MersenneTwister.h"
#include <vector>
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
    ifstream inFile;
    ofstream outFile;
    vector<vector<int>> test(2,vector<int>(2));

    int I = 1;
    int J = 1;

    inFile.open("C.txt", ios::in);
    if (! inFile) {
        cerr << "unable to open file C.txt for reading" << endl;
        return 1;
    }

    for(int i=0; i<I; i++)
        for(int j=0; j<J; j++)
            inFile >> test[i][j];

    outFile.open("results.txt");

    for(int i=0;i<I;i++)
    {
        for(int j=0;j<J;j++)
            outFile<< test[i][j]+1;
        outFile<< endl;
    }

    inFile.close();
    outFile.close();

    return 0;
}



