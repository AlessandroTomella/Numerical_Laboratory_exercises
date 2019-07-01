#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

//Definisco le costanti del problema:

const int Nth = 100000;		//Numero dei lanci
const int BLOCKS = 1000;	//Numero dei blocchi

using namespace std;

int main(){  

//In questa fase inizializzo il generatore di numeri casuali 

 Random rnd;
  
 /*
 Per simulare l'esperimento inziamo prendendo un punto p nel piano in cui la y è scelta in maniera casuale tra 0 e D/2 
 e la x può essere ignorata, come si vedrà, a causa dell'invarainza del problema sotto trasalzioni lungo questo asse. Il punto
 p così trovato identificherà per noi il centro della sbarretta dell'esperimento di Buffon.

 Tirando un secondo punto uniformemente all'interno di un quadrato di lato 2 centrato nell'origine, e accettando solo i punti che cadono
 all'interno della circonferenza goniometrica inscritta possiamo identificare univocamente una direzione nel piano: questo modo di procedere 
 garantisce che la scelta della direzione sia fatta uniformemente tra 0 e 2pigreco. 

 Dall'angolo così identificato possiamo quindi calcolare la coordinata y dell'estremo della sbarretta, a distanza L/2 dal punto 
 p iniziale. L'osservazione importante è che il fatto di aver scelto p come punto medio semplifica moltissimo il codice, come si vedrà, 
 permettendo di concentrarsi effettivamente solo su "un quadrante" del quadrato centrato nell'origine (cioè il quadrante con x e y > 0).
 La scelta di p come punto medio sfrutta al massimo la simmetria del problema per rotazioni di 90 gradi.
 */
 
//Definisco le variabili del problema:

 double L =1;				//Lunghezza della sbarretta
 double D = 2;				//Distanza fra le righe sul pavimento (D > L)
 int Nhit;				//Contatore delle intersezioni
 double y, t, s;			//Variabili per la memorizzazione dei numeri random
 double pi=0, sum=0, sum2=0;		//Variabili per memorizzare il valore di pigreco di ogni esperimento e per la statistica
 ofstream out;				//Variabile di output

 out.open("pi_test.out");
 
 for(int j =0; j < BLOCKS; j++){		//Ciclo sui blocchi
 Nhit =0;					//Riazzero il contatore
 	
	for(int i =0; i < Nth; i++){		//Inizio l'esperimento: ciclo sui lanci della sbarretta
	y = rnd.Rannyu(0, D/2);			//Genero un numero casuale (il centro della sbarretta) tra 0 e D/2 
	t = rnd.Rannyu();			//Genero altri due numeri per identificare una direzione nel piano
	s = rnd.Rannyu();
		while (t*t + s*s > 1){ 		//Ripeto finchè i due punti non stanno DENTRO la circonferenza...
		t = rnd.Rannyu();
		s = rnd.Rannyu();
		}
	y += 0.5*L*sin(atan(s/t));		//...e in tal caso calcolo l'angolo e quindi l'ordinata dell'estremo 
	if (y >= D/2) Nhit++;			//della sbarretta: conto un successo se l'ordinata supera D/2.
	}

 pi = 2*L*Nth/(D*Nhit);				//Dopo l'esperimento faccio la mia stima di pigreco
 sum += pi;					//Aggiorno la media dei blocchi e la varianza
 sum2 += pow(pi, 2);

 //Stampo i risultati in colonna su file e ricomincio.
		
 out << fixed << setprecision(6) << sum/(j+1) << " " << sum2/(j+1) << " " << pow(sum/(j+1),2)<< endl;
 }
 
 out.close();
 return 0;
 }
