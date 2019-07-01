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
const double PI =3.14159265358979323846;//PI greco

using namespace std;

int main(){  

//In questa fase inizializzo il generatore di numeri casuali 

 Random rnd;

//In questa prima fase il calcolo MONTECARLO dell'integrale viene effettuato con un sampling uniforme in [0,1).
//Si procede sempre fornendo varie stime dell'integrale e facendo la statistica sui diversi blocchi.

//Definisco le variabili del problema:

 double x = 0;					//variabile in cui memorizzo il numero casuale
 double mean, m1_mean = 0, m2_mean = 0;		//variabili per l'integrale e i relativi momenti primo e secondo
 double var, m1_var = 0, m2_var = 0;		//variabili per la varianza e i relativi momenti primo e secondo
 ofstream out1, out2;				//variabili di output (uso due file diversi per mean e var)
   
 out1.open("I_unif.out");
 out2.open("V_unif.out");
 

 for (int j =0; j < N; j++){			//Ciclo sul numero di blocchi
 mean=0;					//Riaggiorno a 0 il valore di mean e var nel blocco
 var =0;
	for (int i=0; i<L; i++){		//Ciclo sui punti nel blocco presente
	x = rnd.Rannyu();			//Genero un numero casuale tra 0 e 1
	mean += (PI/2)*cos(PI*x/2)/L;		//Aggiorno i valori di mean e var
	var += pow((PI/2)*cos(PI*x/2) - 1,2)/L;
   	}

 m1_mean += mean;				//Aggiorno i momenti della distribuzione delle mean   
 m2_mean += mean*mean;
 m1_var += var;					//Aggiorno i momenti della distribuzione delle var 
 m2_var += var*var;
 
//Stampo i risultati su file

 out1 << fixed << setprecision(6) << m1_mean/(j+1) << " " << m2_mean/(j+1) << " " << pow(m1_mean/(j+1),2) << endl;
 out2 << fixed << setprecision(6) << m1_var/(j+1) << " " << m2_var/(j+1) << " " << pow(m1_var/(j+1),2) << endl;
 }

 out1.close();
 out2.close();

//Facciamo ora uso dell'IMPORTANCE SAMPLING. Cerchiamo una distribuzione di probabilità che sia "simile" al coseno nell'intervallo [0,1),  
//che sia facile da campionare (i.e. la cui cumulativa sia facilmente invertibile) e che inoltre (visto che tale funzione dovrà 
//apparire al denominatore) non abbia zeri nell'intervallo [0,1), o al più ne abbia agli estremi e in modo da non dare problemi di 
//convergenza. La funzione più semplice che soddisfa questi requisiti è la distribuzione f(x) = a(1-x), con a in modo da garantire la
//normalizzazione. Il procedimento per il resto è uguale:

 double y =0;				//Nuova variabile 
 m1_mean = 0;				//Riazzero i momenti  
 m2_mean = 0;
 m1_var = 0;				 
 m2_var = 0;
 
 out1.open("I_imp.out");		//Apro nuovi file di output
 out2.open("V_imp.out");
 
 for(int j=0; j<N; j++){
 mean=0;					
 var =0;
	for(int i =0; i < L; i++){
	y= rnd.Rannyu();				//Genero un numero uniforme tra 0 e 1
	x = 1 - sqrt(1- y);				//Estraggo quindi un numero dalla mia pdf (la retta)
	
	mean += (PI/4)*cos(PI*x/2)/(1 - x)/L;		//Aggiorno il valore dell'integrale e della varianza
	var += pow((PI/4)*(cos(PI*x/2))/(1 - x) - 1,2)/L;
   	} 
 m1_mean += mean;					//Aggiorno i momenti della distribuzione delle mean   
 m2_mean += mean*mean;
 m1_var += var;						//Aggiorno i momenti della distribuzione delle var 
 m2_var += var*var;
 
//Stampo i risultati su file. 

 out1 << fixed << setprecision(6) << m1_mean/(j+1) << " " << m2_mean/(j+1) << " " << pow(m1_mean/(j+1),2) << endl;
 out2 << fixed << setprecision(8) << m1_var/(j+1) << " " << m2_var/(j+1) << " " << pow(m1_var/(j+1),2) << endl;
 }

 out1.close();
 out2.close();

return 0;
}
