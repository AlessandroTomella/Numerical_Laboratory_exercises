#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

//Definisco le costanti del problema:

const int M = 10000;		//Numero di cammini generati intotale
const int STEPS = 100;		//Numero di passi del singolo cammino
const int BLOCKS = 100;		//Numero di blocchi 
const int L = M/BLOCKS;		//Numero cammini in un blocco

int main(){  

//Inizializzo il generatore di numeri casuali 

 Random rnd;

//Definisco le variabili del problema
 
 double a = 1;
 int c =0, step =0;
 double x[3];
 double r[STEPS + 1], Rmean_block[STEPS + 1], final_mean[STEPS + 1], final_var[STEPS + 1];

//Inizializzo a zero le variabili

 for(int k =0; k < STEPS + 1; k++){
 final_mean[k] = 0;
 final_var[k] = 0;
 }

 ofstream out;
 out.open("data.out");

//In un unico ciclo for genero i cammini, calcolo la distanza dall'origine r, medio all'interno di un blocco e 
//medio infine su tutti i blocchi, calcolando la varianza di r al passo j 

 for(int i =0; i < BLOCKS; i++){					//Ciclo sui blocchi
 for(int k =1; k < STEPS + 1; k++)  Rmean_block[k] = 0;			//Rimetto a 0 la media nel blocco

 	for(int j =0; j < L; j++){					//Ciclo all'interno del blocco
		for(int k =0; k < 3; k++)  x[k] = 0;			//Rimetto a 0 la partenza del cammino e r
		for(int k =0; k < STEPS + 1; k++)  r[k] = 0;

			for(int s =1; s < STEPS + 1; s++){		//Genero un cammino
 			c = rnd.Ranint(0,3);
			step = rnd.Ranint(0,2);
			x[c] += a*step - a*(1-step);
			r[s] = sqrt(pow(x[0],2) + pow(x[1],2) + pow(x[2],2)); //Aggiorno la distanza ad ogni step
			}

		for (int k =1; k < STEPS + 1; k++) Rmean_block[k] += r[k]/L;	//Aggiorno la media dei cammini nel blocco
		}
 	for (int k =1; k < STEPS + 1; k++) {
	final_mean[k] += Rmean_block[k]/BLOCKS;				//Aggiorno le medie dei blocchi
 	final_var[k] += Rmean_block[k]*Rmean_block[k]/BLOCKS;
 	}
 }

//Calcolo la dev standard dalla media di r al passo j

 for (int k =1; k < STEPS + 1; k++) final_var[k] = sqrt(abs(final_var[k] - final_mean[k]*final_mean[k]))/sqrt(BLOCKS);
 
//Stampo in un file i risultati

 for (int k =0; k < STEPS + 1; k++){
 out << fixed << setprecision(6) << final_mean[k] << " " << final_var[k]  << endl;
 }

 out.close();

//In questa seconda fase stampo su file 12 cammini per fare un grafico  


 double walks[12][3], radius = 0;
 out.open("Walks.out");
 
 for(int i =0; i < 12; i++){
	for(int j =0; j < 3; j++) walks[i][j] = 0;
	out << radius <<" ";
 }

 for(int i =1; i < STEPS + 1; i++) {
 out << endl;
	
 	for(int k =0; k < 12; k++){
	c = rnd.Ranint(0,3);
	step = rnd.Ranint(0,2);
	walks[k][c] += a*step - a*(1-step);
	radius = sqrt(pow(walks[k][0],2) + pow(walks[k][1],2) + pow(walks[k][2],2));
	out << fixed << setprecision(6) << radius << " ";
	}
 }

 out.close();
 
//Infine in questa fase salvo un intero cammino su file per potre fare un grafico 3D

 out.open("RW.out");
 for(int k =0; k < 3; k++){  
 	x[k] = 0;
 	out << x[k] << " ";
 }

 for(int j =0; j < 10*STEPS; j++){
	out << endl;
 	c = rnd.Ranint(0,3);
	step = rnd.Ranint(0,2);
	x[c] += a*step - a*(1-step);
	for(int k =0; k < 3; k++) out << x[k] << " ";
 
 }
 
 out.close();		

return 0;

}
