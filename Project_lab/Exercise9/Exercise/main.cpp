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
 int ncities, nstep, init_shuffle;
 double selection_exp, pmutation, pcrossover;
 
 ifstream input;
 input.open("parameters.dat");
 if(input){
 input >> geometry;
 input >> ncities;
 input >> init_shuffle;
 input >> selection_exp;
 input >> pmutation;
 input >> pcrossover;
 input >> nstep;
 }
 else {
 cout << "Could not find file with input parameters. Exit" << endl;
 return 0;
 }
 input.close();


//Creo/Carico la mappa delle città
 TSP Realization(ncities, geometry);

 Realization.LoadCityMap();
 // Realization.GenerateCityMap();
 // Realization.SaveCityMap();

// Parti da una popolazione iniziale pari al quadrato delle città
 int StartingPop = ncities*ncities;

//Dichiaro la variabile Genetic Algorithm 
 GenAlg GA(init_shuffle, selection_exp, pmutation, pcrossover);

//Genero la popolazione iniziale -> 900 individui ciascuno una permutazione di {0, ..., 29}
 vector<Individual> Population;
 for (int i = 0; i < StartingPop; i++){
// Genera un individuo
 Individual ind(Realization, ncities, geometry);
// Disordinane gli elementi 
 GA.Shuffle(ind);
// Controlla la validità dell'individuo
 if(ind.SelfCheck()) Population.push_back(ind);

 else {
 	cout << "Error in generating starting sequence /////////////////////////////////" << endl;
 	return 0;
 	}
  
 }
 
 double aveFitness;
 ofstream output;
 if (geometry) output.open("WorkInProg/Circular/Fitness.dat");
 else output.open("WorkInProg/Square/Fitness.dat");
 
// Inizia l'evoluzione della popolazione
 for(int k = 1; k <= nstep; k++){

 cout << "//////////////////////// Generazione numero : " << k << endl;

 GA.Hierarchy(Population);		// Ordina la popolazione attuale
 Population.at(0).PrintInd(k); 		// Salva il miglior individuo

 aveFitness = 0;
 for(int i = 0; i < Population.size(); i++){	// Calcola l'average fitness della popolazione
	aveFitness += Population.at(i).GetFitness()/(double)Population.size();
 }
 
 output << aveFitness << endl;		// Stampa l'average fitness

 Population = GA.Evolve(Population);	// Chiama la funzione principale del codice, che gestisce un passo evolutivo

 for (int i = 0; i < Population.size(); i++){	// Controlla la popolazione
	
	if(!Population.at(i).SelfCheck()){
		cout << "Error in mutated population ////////////////////////////" << endl;
		return 0;		
	}
 }

}

// Ordina la popolazione finalee mostra l'individuo con la sua fitness
GA.Hierarchy(Population);
cout << "** Miglior percorso dell'ultima generazione **" << endl;
Population.at(0).ShowInd();
cout << "Lunghezza: "<< Population.at(0).GetFitness() << endl;


return 0;
}
