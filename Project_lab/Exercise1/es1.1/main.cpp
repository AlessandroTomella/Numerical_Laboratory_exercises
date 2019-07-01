#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

//Definisco le costanti del problema:

const int M=100000;		//Numero totale punti generati
const int N=100;		//Numero dei blocchi
const int L=M/N;		//Numero dei punti in un blocco
const int INT=100;		//Numero intervalli in cui suddivido [0,1)
const int CHI = 10000;		//Numero di punti generati per il test del chi_quadro

using namespace std;

int main(){  

 Random rnd;
  
//Nella prima parte dell'esercizio stimo il valore dei due integrali con un metodo MONTECARLO, generando L punti casuali uniformemente
//tra 0 e 1 e calcolando il valor medio della funzione integranda su tali punti. La procedura viene ripetuta N volte e i risultati 
//vengono utilizzati per fare la statistica sui valori ottenuti: vengono calcolati i momenti primo e secondo in funzione 
//del numero di simulazioni effettuate fino ad un determinato punto. Il valore finale (cioè quello ottenuto dopo tutte le N simulazioni) 
//sarà la nostra migliore stima del'integrale e dell'errore su questa quantità.

//Definisco le variabili del problema. 
//COMMENTO: Con i termini "mean" e "var" mi riferisco sempre al valore 
//          del primo e del secondo integrale dell'esercizio! 

 double x = 0;					//variabile in cui memorizzo il numero casuale
 double mean, m1_mean = 0, m2_mean = 0;		//variabili per la stima del primo integrale (mean)
 double var, m1_var = 0, m2_var = 0;		//variabili per la stima del secondo integrale (var)
 ofstream out1, out2;				//variabili di output (uso due file diversi per mean e var)
 
 out1.open("mean.out");
 out2.open("variance.out");
	
//Ciclo che produce i risultati

 for (int j =0; j < N; j++){			//Ciclo sul numero di blocchi
 mean=0;					//Riaggiorno a 0 il valore di mean e var nel blocco
 var =0;
	for (int i=0; i<L; i++){		//Ciclo sui punti nel blocco presente
	x = rnd.Rannyu();			//Genero un numero casuale tra 0 e 1
	mean += x/L;				//Aggiorno i valori di mean e var
	var += (x - 0.5)*(x - 0.5)/L;
	}
 m1_mean += mean;				//Aggiorno i momenti della distribuzione delle mean   
 m2_mean += mean*mean;
 m1_var += var;					//Aggiorno i momenti della distribuzione delle var 
 m2_var += var*var;
	
//Stampo in colonna i valori momento primo, del momento secondo e del momento primo al quadrato al passo attuale
//cioè tenendo conto di tutti i blocchi finora generati.

 out1 << fixed << setprecision(6) << m1_mean/(j+1) << " " << m2_mean/(j+1) << " " << pow(m1_mean/(j+1),2) << endl;
 out2 << fixed << setprecision(6) << m1_var/(j+1) << " " << m2_var/(j+1) << " " << pow(m1_var/(j+1),2) << endl;
 }
  
 out1.close();
 out2.close();

//Per non complicare le cose il test del chi-quadro viene effettuato in un secondo momento, e quindi inevitabilmente (a meno di
//ri-inizializzare il generatore) con numeri random diversi da quelli usati per valutare i due integrali. 
//La procedura qui sarà quella di suddividere l'intervallo [0,1) in INT intervalli uguali, generare CHI numeri casuali nell'intervallo 
//[0,1), contare quanti punti cadono in ogni intervallino e con questi valori fare una prima stima del parametro chi_quadro. Anche qui
//ripetendo la procedura N volte si potrà fare la statistica sui risultati.
 
//Anche qui definisco le mie variabili.

 double chi, m1 = 0, m2 = 0;			//Variabili in cui memorizzo il valore chi-quadro e delle medie cumulative
 int n_points[N];				//Variabile che conta quanti punti cadono nell'intervallo i-esimo

 out1.open("chi_values.out");

 for(int j =0; j < N; j++){			//Ciclo sul numero di blocchi
 	chi =0;					//Riaggiorno a zero il valore di chi e del vettore dei punti contati
 	for(int r=0; r < INT; r++) n_points[r] = 0; 
	
	for(int i =0; i < CHI; i++){		//Ciclo sui punti nel blocco presente
		x = rnd.Rannyu();			//Genero un numero casuale tra 0 e 1
		for(int k=0; k < INT; k++){	//Ciclo sugli intervallini e segno in quale cade il punto x generato
			if (x>=(double)k/(double)100 && x<(double)(k+1)/(double)100) n_points[k]++;
		}
	}
  	for(int s=0; s < INT; s++){		//Calcolo il valore del chi-quadro in questo blocco
		chi += pow(n_points[s] - (float)CHI/INT,2);
 	}
 chi /= (double)CHI/(double)INT;
 
 m1 += chi;					//Aggiorno i momenti della distrubuzione dei chi
 m2 += chi*chi;				
 
//Stampo i risultati su file

 out1 << fixed << setprecision(6) << m1/(j+1)  << " " << m2/(j+1) << " " << pow(m1/(j+1),2) <<endl;
 }
 
 out1.close();
 
return 0;
}
