/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "funz.h"

using namespace std;
 
Random :: Random(){
 int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
 
 ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
return;
}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

/*************************************************************************************/
/*************************************************************************************/

/*     Metodi classe RW    */

/*************************************************************************************/
/*************************************************************************************/

RW :: RW (){

count = 0;	// Mosse accettate
EQ_STEPS = 50;  // Passi random walk per equilibrare 

// Il passo del random walk è stato scelto in modo da avere l'accettazione Metropolis al 50%.
// Questo è stato fatto solo per i due orbitali richiesti dall'esercizio, per gli altri si è usato
// un valore generico.
if ( orb.get_n() == 1 && orb.get_l() == 0 && orb.get_m() == 0) a = 1.2;
else if ( orb.get_n() == 2 && orb.get_l() == 1 && orb.get_m() == 0) a = 2.6;
else a = 3;


cout << "Choose a starting point for sampling." << endl;
cout << "x = " ;
cin >> x[0];
cout<< endl;
cout << "y = " ;
cin >> x[1];
cout<< endl;
cout << "z = " ;
cin >> x[2];
cout<< endl;

}

RW :: ~RW(){}

// Equilibrazione, esegue alcuni passi del random walk a partire dalla posizione attuale
void RW :: Equilibrate(){

char a;
PrintPoint();
cout << "Do you want to integrate a few steps before starting the sampling? [y/n]...";
cin >> a;
while(a != 'y' && a!= 'n'){
cout << "Invalid answer! Must be either 'y' or 'n', try again..."<< endl;
cin >> a;
}
while(a == 'y'){
	for(int i = 0; i < EQ_STEPS; i++) U_Step();
	PrintPoint();
	cout << "Do you want to integrate a few more steps before starting the sampling? [y/n]...";
	cin >> a;
}

count = 0;
return ;

}

// STampa la posizione attuale a video
void RW :: PrintPoint(void){

cout << "Present position: " << endl;
cout << " x = " << x[0] << endl;
cout << " y = " << x[1] << endl;
cout << " z = " << x[2] << endl;  
}

//Restituisce le coordinate del punto attuale
double RW :: GetCoord(int n){
while( n != 0 && n != 1 && n != 2){
cout << "Index out of bound for 3-dimensional vector! Please enter tthe index for the wanted coordinate [0 for x, 1 for y, 2 for z]..." << endl;
cin >> n;
cout << endl;}

return x[n];
}

// Esegue il passo del random walk nel continuo 
// con direzione scelta uniformemente nell'angolo solido
void RW :: U_Step(void){

double theta, phi, p, pnew;
double y[3];
theta = acos(1-2*rnd.Rannyu());		// Angolo zenitale
phi = rnd.Rannyu(0, 2*PI);		// Angolo azimutale

p = orb.Eval(x);			// Valuto la probabilità della configurazione attuale				

y[0] = x[0] + a*sin(theta)*cos(phi);	// Salvo la nuova posizione
y[1] = x[1] + a*sin(theta)*sin(phi);	// nella variabile 
y[2] = x[2] + a*cos(theta);		// temporanea y
pnew = orb.Eval(y);			// e valuto la probabilità associata a questa nuova configurazione

if (Accept(pnew, p)){			// Chiamo la funzione di accettazione
 x[0] = y[0];				// In caso affermativo sovrascrivo x
 x[1] = y[1];
 x[2] = y[2];
 count++;				// e aggiorno il contatore delle mosse accettate
 }

}

// Simile alla funzione precedente, ma il passo è scelto
// estraendo da una gaussiana centrata in 0 e con varianza passo_uniforme/3
void RW :: G_Step(void){

double x_step, y_step, z_step, p, pnew;
double y[3];
x_step = rnd.Gauss(0, a/sqrt(3));
y_step = rnd.Gauss(0, a/sqrt(3));
z_step = rnd.Gauss(0, a/sqrt(3));

p = orb.Eval(x);				

y[0] = x[0] + x_step;
y[1] = x[1] + y_step;
y[2] = x[2] + z_step;

pnew = orb.Eval(y);

if (Accept(pnew, p)){
 x[0] = y[0];
 x[1] = y[1];
 x[2] = y[2];
 count++;
 }

}


bool RW :: Accept(double x, double y){
if((x > y)||(rnd.Rannyu() < x/y)) return 1;
else return 0;
}


int RW :: GetCount(void){
return count;
}

// Calcolae restituisce la distanza dall'origine della posizione corrente 
double RW :: GetRadius(void){
 double r = 0;
 for(int i = 0; i < 3; i++) r += x[i]*x[i];
 r = sqrt(r);

return r;
}

// Numero quantico pricnipale
int RW :: Getn(){
return orb.get_n();
}

// Numero quantico azimutale
int RW :: Getl(){
return orb.get_l();
}

// Numero quantico magnetico
int RW :: Getm(){
return orb.get_m();
}


/*************************************************************************************/
/*************************************************************************************/

/*     Metodi classe Orbital   */

/*************************************************************************************/
/*************************************************************************************/

// Nel costruttore viene chiesto all'utente di inserire i numeri quantici dell' orbitale
// che vuole campionare, con alcuni controlli su tale input
Orbital :: Orbital (){

cout << "Enter the quantum numbers that define the orbital." << endl;
cout << "Principal quantum number: n = " ;
cin >> n ;
cout<< endl;
while(n < 1){
	cout << "Invalid choice for quantum number n: n must be a positive integer! Choose a different n." << endl;
        cout << "n = " ;
	cin >> n;
	cout << endl;	
	}
cout << "Azimuthal quantum number: for computational reasons only 3 values for l will be accepted, namely 0, 1 and 2." << endl;
cout << "Check also that your choice for l satisfies the condition n-l-1 >= 0. " << endl;
cout << " l = " ;
cin >> l;
while((n-l-1 < 0) || (l != 0 && l != 1 && l != 2)) {
	cout << "Invalid choice for quantum number l, be sure that:" << endl;
	cout << "- l is either 0, 1 or 2 " << endl;
	cout << "- n-l-1 is a non-negative quantity" << endl;
	cout << "l = ";
	cin >> l;
	cout << endl;
	}
if(l == 0) {
	m = 0;
	cout << "Magnetic quantum number: m = " << m << endl;
	return;
	}
else if (l == 1){
	cout << "Magnetic quantum number: only 3 values for m will be accepted, namely 0 and +/- 1." << endl;
	cout << " m = " ;
	cin >> m;
	cout << endl;
	while(m != 0 && m*m != 1) {
		cout << "Invalid choice for quantum number m: m must be either 0 or +/- 1! Choose a different m." << endl;
		cout << "m = ";
		cin >> m;
		cout << endl;
		}
	return;
	}
else {
	cout << "Magnetic quantum number: only 5 values for m will be accepted, namely 0, +/- 1 and +/- 2." << endl;
	cout << " m = " ;
	cin >> m;
	cout<< endl;
	while(m != 0 && m*m != 1 && m*m != 4) {
		cout << "Invalid choice for quantum number m: m must be either 0, +/- 1or +/- 2! Choose a different m." << endl;
		cout << "m = ";
		cin >> m;
		cout << endl;
		}
	return;
	}
return;
}

Orbital :: ~Orbital(){}

int Orbital :: get_n(){
return n;
}

int Orbital :: get_l(){
return l;
}

int Orbital :: get_m(){
return m;
}

// Calcolo del fattoriale
int Orbital :: Fact(int j){
 int f = 1;
 for(int i = 1; i < j+1; i++) f *= i;
 return f;
}

//Calcolo del coefficiente binomiale
double Orbital :: Bin(int j, int k){
 if(k > j) {
 cout << "Error in binomial function" << endl;
 return 0;}

 return (double)Fact(j)/(double)(Fact(k)*Fact(j-k));
}

// Fattore di normalizzazione della funzione d'onda
double Orbital :: Norm(void){
 int j = n - l - 1;

 return pow(2.0/(double)n, 3) * Fact(j)/(double)(2.0 * n * Fact(n + l));
}

// Restituisce il polinomio di Laguerre valutato in x
double Orbital :: Laguerre(double x){
 int j = n - l - 1;
 int a = 2*l + 1;
 double y = 0;
 if(m < 0) {
 cout << "Error in Laguerre polynomial" << endl;
 return 0;
 }
 else {
	for(int i = 0; i < j+1; i++){
	y += pow((-1), i) * Bin(j + a, j - i) * pow(x, i) / (double)Fact(i); 
	}
 }

 return y*y;
}

// Restituisce il polinomio di Legendre (armonica sferica) valutato in x ( = cos(theta)) 
double Orbital :: Harmonic(double x){

 if(m == 0){
	if(l == 0) return 0.25/PI;	
	if(l == 1) return 0.25 * 3/PI * x * x;		
	if(l == 2) return 0.25 * 0.25 * 5 / PI * pow((3 * x * x - 1), 2);
 }
 if(m == 1 || m == -1){
	if(l == 1) return 0.25 * 3 /(2*PI) * (1 - x * x);	
	if(l == 2) return 0.25 * 15 /(2*PI) * x * x * ( 1 - x * x);	
 }
 if(m == 2 || m == -2){
	if(l == 2) return 0.25 * 0.25 * 15 / (2*PI) * pow(1 - x * x, 2);	
 }
else{ cout << "Error in Harmonic funtion!" << endl; }
return 0;
}

// Calcola il modulo quadro della funzione d'onda 
double Orbital :: Eval(double r[3]){
 double R = 0;
 for(int i = 0; i < 3; i++) R += r[i]*r[i];
 R = sqrt(R);
 double cos_theta = r[2]/R;

 return Norm() * Laguerre(2*R/(double)n) * Harmonic(cos_theta) *  pow((2*R/(double)n), 2*l) * exp(- 2*R/(double)n);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
