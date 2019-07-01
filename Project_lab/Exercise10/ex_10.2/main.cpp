#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <bits/stdc++.h> 
#include "vector"
#include "funz.h"
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]){
 
MPI::Init(argc, argv);
	
 int world_size = MPI::COMM_WORLD.Get_size();
 int world_rank = MPI::COMM_WORLD.Get_rank();
 
// Questa parte è un simulated annealing uguale in tutti processori

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

 SimAnneal SA(init_shuffle, pmutation, world_rank);

//Genero l'individuo iniziale
 Individual ind(Realization, ncities, geometry);

 SA.Shuffle(ind);
 ind.PrintInd(0);
 
 int count = 0;
 while(beta <= 1000){
	
	for(int k = 1; k <= nstep; k++){
	SA.Metropolis(ind, beta);
	}
 
 beta += beta_incr; 
 SA.ResetAccepted();
 }

// Parte parallelizzata per raccogliere i risultati

 	double irecv[4];
	for(int i=0;i<4;i++)
		irecv[i]=0;
	
	double isend = ind.CostFunc();;
	cout << "Fitness del processo " << world_rank << " = " << isend << endl;
	MPI_Gather(&isend,1,MPI_DOUBLE,irecv,1,MPI_DOUBLE,0,MPI::COMM_WORLD);
	int who;
	int *ID = new int[ncities];
	if(world_rank==0){
		cout<<"Gathered Fitnesses: " << irecv[0] << " " << irecv[1]<<" "<<irecv[2]<<" "<<irecv[3]<<endl;	
	
	double best; 
	best = irecv[0];
	for(int i =1; i < 4; i++){
		if  (irecv[i] <= best){
		best = irecv[i];
		who = i;
		}
	}
	cout << "Best: " << best << " in process " << who << endl;
}

MPI_Barrier(MPI_COMM_WORLD);
MPI_Bcast(&who, 1, MPI_INT, 0, MPI_COMM_WORLD);
if(who == world_rank) {
	cout << "It's me! Process " << world_rank << endl;
	ind.ShowInd();
	ind.PrintInd(1);	 	
}
MPI::Finalize();

return 0;
}
