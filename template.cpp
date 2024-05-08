#include <iostream>             // les 4 1eres lignes c'est pour appeler des librairies pour sauvegarder des fichiers et utiliser des fonctions mathematiques
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include "MersenneTwister.h"    // les deux c'est pour mettre un commentaire cette ligne te permet d'utiliser mersenne

using namespace std;

int main()           // ça c'est ton programme principal ça commence tout le temps comme ça
{   
    ofstream files ("<k2>-<k>.dat");
    files << 0 << " " << 0 << endl;
    files.close();

    ofstream fil ("<r2>+3<r>2.dat");
    fil << 0 << " " << 0 << endl;
    fil.close();

    ofstream  fyl ("<r2>+3<r>2-<r>.dat");
    fyl << 0 << " " << 0 << " " << 1 << " " << 1 << endl;
    fyl.close();


    
    MTRand mtrand1;          // ça c'est pour pouvoir utiliser twister il faut le déclarer

    int const N(10000000);  // tu déclare toutes tes variables en leur assignant ou non une valeur les () ou = c'est pareil
    int N1;
    double t(100);
    double ta(0);
    double tau;
    
    double alpha(0.7);    // utilise des entiers si tu peux c'est plus rapide à comparer et à faire des opérations dessus
    double dmax(0);
    double r(0.5);
    double r_avg(0);
    double deltaT = 100;
    double sum(0);
    double Dmax(0);
    double c_0(1);
    double betha(3);
    double eps = 0.001;
    double gamma(1);
    double r2(0);
    double A;
    double c_avg;

    double *act = new double[N];      // ça c'est pour faire un vecteur de taille N en demandant de la mémoire à l'ordi, en fait c'est un pointeur mais ça s'utilise comme un vecteur mais il faut bien initialiser tes valeurs car la case mémoire contient éventuellement déjà qq chose 
    int *count = new int[N];
   
    
    for (int i = 0; i < N; i++)         // une boucle for
    {
        act[i] = 0;
        count[i] = 0;
        
        while (tau < ta + v*deltaT)  // dans les conditions si tu veux mettre un égal il faut écrire ==, après tu as <= (inférieur ou égal) >= et != (différent de)
        {
            if (tau > ta)
            {
                                                 
            }
        }
    }
    
    // fais gaffe si tu divises par un entier il te rend la division euclidienne donc faut écrire
    
    1/(double)N // par exemple
           

  // commenter un bloc          
        
 /*       
        for (int i = 0 ; i < N ; i++)
        {
            r2 += pow(count[i]-r_avg,2)/(double)N;
        }

 */     

            
	int n = mtrand1.randInt(N-1);   // entier aléatoire entre 0 et N-1 inclus    
	r = mtrand1.randDblExc(); // reel aleatoire entre 0 et 1 exclus  c est expliqué dans le fichier mersenne.h
        
    
      
    ofstream pop("degDistr.dat");        // ça c'est pour enregistrer dans un fichier txt, tu peux utiliser ce que tu veux à la place de pop

    for (int compteur = 0 ; compteur < D; compteur++)
    {
        pop << compteur << " " << act[compteur] << " " << " " << count[compteur] << endl;
    }
    pop.close();

// le endl; c'est pour sauter une ligne donc dans ton fichier t'auras

/*

 0  act[0]  count[0]
 1  act[1]  count[1] 
 ....


*/    
    
    delete [] count;    // pour rendre la mémoire que tu as empruntée
    delete [] act; 
    cout << r_avg << endl;    // pour afficher la valeur de r_avg à l'écran
}



