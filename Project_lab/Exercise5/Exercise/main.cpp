#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include "vector"
#include "funz.h"

//constants
const int M = 1000000;		//Numero passi del random walk

const int BLOCKS = 100;		//Numero di blocchi 
const int L = M/BLOCKS;		//Numero cammini in un blocco
const int PRINT = 100;		//Numero passi tra successiva stampe su file

using namespace std;

void UpdateTable(vector<vector<int>> &v, vector<int> last);	// Tabella che memorizza gli orbitali già campionati
void DisplayTable(vector<vector<int>> );			// e la stampa a video a ogni iterazione


int main(){  
 
 bool sample_again = 1;			// Controlla uscita dal programma
 int conf = 1;				// Numero di sampling eseguiti
 
 double R, R_ave, R_var;		// Distanza dall'origine e valori medi
 int N;					//  	
 
 vector<vector<int>> table;		// Tabella con orbitali campionati
 
 ofstream out;				// Variabili di stampa su file
 ofstream out_ave;

while(sample_again){

 cout << "************* STARTING THE SAMPLING NUMBER " << conf << " ****************" << endl << endl;
 
 R = 0;			// Riazzero tutte le variabili per iniziare 
 R_ave = 0;		// un nuovo campionamento
 R_var = 0;
 N = 0;
 
 out.open("WorkInProg/State_2p.out");	// Salvo punti campionati in 3D
 out_ave.open("WorkInProg/R_2p.out");	// Salvo valori medi nei blocchi della distanza dall'origine
 
 if(conf > 1) DisplayTable(table);			// Mostra la tabella di ciò che è già stato fatto
 
 RW walker;				// Nel costruttore di walker vengono inizializzate tutte le variabili
					// e viene scelto un orbitale da campionare
 
 out << walker.Getn() << " " << walker.Getl() << " " << walker.Getm() << " " << endl;  // Stampo nella prima riga del file i n° quantici  
 vector<int> numb = { {walker.Getn(), walker.Getl(), walker.Getm()} }; 		// Aggiorno la tabella con il nuovo orbitale

 walker.Equilibrate();				// Fase di equilibrazione (alcuni passi del random walk)
 cout << "Starting the sampling at the ";
 walker.PrintPoint();				// Stampa coordinate del punto attuale

 for(int i = 0; i < M; i++){		
 
 	if(i%PRINT == 0){			// Stampa su file le coordinate attuali
		out << fixed << setprecision(6) << walker.GetCoord(0) << " " << walker.GetCoord(1)  << " " << walker.GetCoord(2)  << endl;
 	}			

 	R += walker.GetRadius()/(double)L;	// Aggiorna la media della distanza dall'origine nel blocco 
 
	if((i+1)%L == 0){			// Aggiorno le medie finali nel blocco	
		N++;
		R_ave += R;
		R_var += R*R;
	
		// Salvo su file 3 colonne corrispondenti a: R medio, R quadro medio, R medio quadro
		out_ave << R_ave/(double)N << " " << R_var/(double)N << " " <<  pow(R_ave/(double)N, 2) << endl;

		R = 0;				// Riazzero la media nel blocco
		}

 	//walker.U_Step();		// Mossa del random walk, con accettazione secondo Metropolis
	walker.G_Step();		
 }
			
 cout << "Accepted to proposed move ratio: " << (double)walker.GetCount()/(double)M << endl;
 
 out.close();		// Chiudo l'ofstream
 out_ave.close();
 
 conf++;			// Aggiorno il numero di orbitali campionati
 UpdateTable(table, numb);	// Aggiornio la tabella degli orbitali campionati


 cout << "Do you want to sample another orbital? The old results will not be deleted!"<< endl;
 cout << "Type any number if you want to sample again, 0 if you want to exit the program." << endl;
 cin >> sample_again;
}
 
return 0;
}

// Aggiunge una riga alla fine della tabella
void UpdateTable (vector<vector<int>> &v, vector<int> last){

v.push_back(last);
return;
}


// Stampa a video la tabella
void DisplayTable(vector<vector<int>> grid){

cout << "Already sampled orbitals:" << endl;
cout << "      " << "n  " << "l  " << "m  " << endl; 
cout << "     -----------" << endl;
for (const vector<int> v : grid )
{
   cout << "      " ;
   for ( int x : v ) cout << x << "  ";
   cout << endl;
}
return;
}


