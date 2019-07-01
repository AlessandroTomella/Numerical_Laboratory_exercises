#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

//Definisco le costanti del problema:

const int M = 10000;				//Numero punti nell'istogramma
const double pi = 3.14159265358979323846;	//PI greco

using namespace std;

int main(){  

//In questa fase inizializzo il generatore di numeri casuali 

 Random rnd;
  
//In questo codice si vogliono generare vari sample di random variabili distribuite secondo 3 diverse distribuzioni di probabilità
//e studiare inoltre come sono distribuite le random variabili Medie, ottenute come medie su N valori delle random variabili iniziali,
//con N che assume 4 i valori 1, 2, 10, 100.
 
//Definisco le variabili del problema

 ofstream out1, out2, out3;		//Variabili di output, 3 diversi file per le 3 distribuzioni
 double Sum_U[4], Sum_E[4], Sum_L[4];	//3 Vettori (1 per ogni pdf) di 4 elementi ciascuno, corrispondenti a 4 diverse variabili Media
 double x =0, exp = 0, lor = 0;		//Variabili in cui memorizzo il numero generato secondo la distribuzione di probabilità
 int N[4] = {1, 2, 10, 100};		//Vettore che determina le 4 variabili Media
 
 out1.open("ST_dice.out");
 out2.open("Exp_dice.out");
 out3.open("Lor_dice.out");

 for(int k =0; k < M; k++){ 				//Ciclo sul numero di punti che voglio nell'istogramma finale
	for (int j =0; j < 4; j++){ 			//Ciclo sull'estremo N delle variabili Medie
	Sum_U[j] = 0;					//Riaggiorno a 0 il valore delle variabili Media
	Sum_E[j] = 0;
	Sum_L[j] = 0;
		for(int i = 0; i < N[j]; i++){		//Sommo N valori, dove N dipende dal ciclo
		x = rnd.Rannyu();			//Estraggo dalla distribuzione uniforme
		while (x == 0) {x = rnd.Rannyu();}  	//Assicuro che x sia diverso da 0, per evitare errori nel calcolo 
							//della lorentziana
	        exp = -log(1-x); 			//Estraggo da un esponenziale
		lor = tan((x - 1/2)*pi);		//Estraggo da una lorentziana
		
		Sum_U[j] += x;				//Aggiorno il valore della variabile Media di N
		Sum_E[j] += exp;
		Sum_L[j] += lor;
		}
 	}

//Stampo i risultati su 3 file diversi in modo che su ciascuno appaiano 4 colonne, corrispondenti alle 4 variabili Media per diverso N,
//ciascuna costituita da M valori diversi, che saranno quelli da plottare nell'istogramma.

 out1 << fixed << setprecision(6) << Sum_U[0]/N[0] << " " << Sum_U[1]/N[1] << " "<< Sum_U[2]/N[2] << " "<< Sum_U[3]/N[3] <<endl;
 out2 << fixed << setprecision(6) << Sum_E[0]/N[0] << " " << Sum_E[1]/N[1] << " "<< Sum_E[2]/N[2] << " "<< Sum_E[3]/N[3] <<endl;
 out3 << fixed << setprecision(6) << Sum_L[0]/N[0] << " " << Sum_L[1]/N[1] << " "<< Sum_L[2]/N[2] << " "<< Sum_L[3]/N[3] <<endl;
 }
	 
 out1.close();
 out2.close();
 out3.close();
 
 return 0;

}



