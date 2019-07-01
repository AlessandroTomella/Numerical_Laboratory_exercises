#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

//Definisco le costanti:

const int M=100000;			//Numero totale punti generati
const int N=100;			//Numero dei blocchi
const int L=M/N;			//Numero dei punti in un blocco

const double S_0 = 100;			//Valore di mercato dell'asset al tempo iniziale
const double T = 1;			//Data di expiry dell'opzione
const double K = 100;			//Strike price
const double rfir = 0.1;		//Risk-free interest rate
const double vol = 0.25;		//Volatility
const double t_steps = 100;		//Numero di intervalli temporali in cui discretizzare il processo

double maxp(double x, double y){	//Funzione che restituisce il massimo tra due valori
 if (x>=y) return x;
 else return y;
};

using namespace std;

int main(){  

//In questa fase inizializzo il generatore di numeri casuali 

 Random rnd;


//Come prima cosa calcoliamo il prezzo dell'opzione simulando l'andamento dell'asset price direttamente al tempo finale T.

//Definisco le variabili del problema:

 double S, call_price, put_price;			//Variabili in cui salvo il prezzo dell'asset, dell'opzione call e dell'opzione put
 double c_mean = 0, c_var = 0, p_mean = 0, p_var = 0;	//Variabili per la statistica nei blocchi
 ofstream out1, out2;				 		//Variabili di output
   
 out1.open("Dir_Call.out");
 out2.open("Dir_Put.out");
 
 for (int j =0; j < N; j++){							//Ciclo sul numero di blocchi
 call_price = 0;
 put_price = 0;
	for (int i=0; i<L; i++){						//Ciclo sui punti nel blocco presente
	S = S_0*exp((rfir - 0.5*vol*vol)*T + vol*rnd.Gauss(0, sqrt(T)));	//Genero il prezzo dell'asset al tempo T
	call_price += exp(-rfir*T)*maxp(0, S - K)/L;				//Ricavo i relativi prezzi di un'opzione call
	put_price += exp(-rfir*T)*maxp(0, K - S)/L;				//e di un'opzione put
   	}
 c_mean += call_price;								//Aggiorno i momenti della distribuzione delle call   
 c_var += call_price*call_price;
 p_mean += put_price;								//Aggiorno i momenti della distribuzione delle put   
 p_var += put_price*put_price;

//Stampo i risultati su file

 out1 << fixed << setprecision(6) << c_mean/(j+1) << " " << c_var/(j+1) << " " << pow(c_mean/(j+1),2) << endl;
 out2 << fixed << setprecision(6) << p_mean/(j+1) << " " << p_var/(j+1) << " " << pow(p_mean/(j+1),2) << endl;
 }

 out1.close();
 out2.close();

//Vogliamo ora ripetere il calcolo, discretizzando l'intervallo temporale T in t_steps intervalli intermedi, e simulando
//l'andamento dell'asset price nel tempo: poichè la discretizzazione è esatta il risultato finale del calcolo deve essere
//lo stesso ottenuto con il metodo "diretto" precedente. 

 c_mean = 0;				//Riaggiorno a zero le variabili								
 c_var = 0;
 p_mean = 0;								   
 p_var = 0;

 out1.open("Disc_Call.out");
 out2.open("Disc_Put.out");
 
 for (int j =0; j < N; j++){						//Ciclo sul numero di blocchi
 call_price = 0;
 put_price = 0;
	for (int i=0; i<L; i++){					//Ciclo sui punti nel blocco presente
	S = S_0;
		for(int k = 0; k < t_steps; k++){			//Simulo l'andamento dell'asset price nel tempo	
		S = S*exp((rfir - 0.5*vol*vol)*(T/t_steps) + vol*rnd.Gauss(0, 1)*sqrt(T/t_steps));	 
		}
	call_price += exp(-rfir*T)*maxp(0, S - K)/L;			//Ricavo i relativi prezzi di un'opzione call
	put_price += exp(-rfir*T)*maxp(0, K - S)/L;			//e di un'opzione put
   	}
 c_mean += call_price;							//Aggiorno i momenti della distribuzione delle call   
 c_var += call_price*call_price;
 p_mean += put_price;							//Aggiorno i momenti della distribuzione delle put   
 p_var += put_price*put_price;

//Stampo i risultati su file

 out1 << fixed << setprecision(6) << c_mean/(j+1) << " " << c_var/(j+1) << " " << pow(c_mean/(j+1),2) << endl;
 out2 << fixed << setprecision(6) << p_mean/(j+1) << " " << p_var/(j+1) << " " << pow(p_mean/(j+1),2) << endl;
 }

 out1.close();
 out2.close();


return 0;
}
