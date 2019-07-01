#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <bits/stdc++.h> 
#include "vector"
#include "funz.h"

using namespace std;

int main(){
 
 bool geometry;
 int ncities, init_shuffle, nstep;
 double pmutation, beta, beta_incr;
 
 ifstream input;
 input.open("parameters.dat");
 if(input){
 input >> geometry;
 input >> ncities;
 input >> init_shuffle;
 input >> pmutation;
 input >> beta;
 input >> beta_incr;
 input >> nstep;
 }
 else {
 cout << "Could not find file with input parameters. Exit" << endl;
 return 0;
 }
 input.close();

 Random rnd;
//Creo/Carico la mappa delle città
 TSP Realization(ncities, geometry);

 //Realization.LoadCityMap();

 Realization.GenerateCityMap(rnd);
 Realization.SaveCityMap();

// Creo l'oggetto Simulated Annealing
 SimAnneal SA(init_shuffle, pmutation);

//Genero l'individuo iniziale
 Individual ind(Realization, ncities, geometry);
 SA.Shuffle(ind);
 ind.PrintInd(0);

 ofstream output;
 if (geometry) output.open("Circle/Annealing.dat");
 else output.open("Square/Annealing.dat");
 
 int count = 0;
 while(beta <= 1000){
	cout << "////////////////////////////////////////" << endl;
	cout << "Simulating a temperature of " << 1.0/beta << " for " << nstep << " MC steps..." << endl;
	for(int k = 1; k <= nstep; k++){
	SA.Metropolis(ind, beta);
	}
 //ind.ShowInd();
 cout << "Acceptace ratio " << SA.GetAccepted()/(double)nstep << endl << endl;
 if(SA.GetAccepted() > 0) {
 	count++;
 	ind.PrintInd(count);	// salvo indivduo su file solo se è diverso dal precedente
 }
 
 output << beta << " " << ind.CostFunc() << endl;
 beta += beta_incr; 
 SA.ResetAccepted();
 }
 
 output.close();
 

return 0;
}
