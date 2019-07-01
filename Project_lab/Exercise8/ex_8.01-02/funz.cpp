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
  // string property;
   if (input.is_open()){
     // while ( !input.eof() ){
       //  input >> property;
         //if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            SetRandom(seed,p1,p2);
        // }
     // }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
return;
}

Random :: ~Random(){
cout << "Saving the Seed..." << endl;
SaveSeed();
}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.in");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
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

 count = 0;
 a = 5;
 SetPoint(1.0);
/* cout << "Choose a starting point for sampling." << endl;
 cout << "x = " ;
 cin >> x;
 while(x >=20 || x <= -20){
	cout << "Invalid input, it's too big!(That's what she said)" << endl;
	cout << "Choose again:"<< endl;
	cin >> x;
 }
 cout << endl;
/*/
}

RW :: ~RW(){}

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
	for(int i = 0; i < EQ_STEPS; i++) Mrt_Step();
	PrintPoint();
	cout << "Do you want to integrate a few more steps before starting the sampling? [y/n]...";
	cin >> a;
}
return ;

}


void RW :: SetPoint(double r){

x = r;
return;
}

void RW :: PrintPoint(void){

cout << "Present position: " << endl;
cout << " x = " << x << endl;

}

double RW :: GetPoint(void){

return x;
}

void RW :: ResetCount(){

count = 0;
}


void RW :: Mrt_Step(void){

double y, p, pnew;

p = wf.Eval(x)*wf.Eval(x);				

y = x + a*(rnd.Rannyu()-0.5);

pnew = wf.Eval(y)*wf.Eval(y);

if (Accept(pnew, p)){
 x=y;
 count++;
 }

}
/*
void RW :: G_Step(void){

double x_step, y_step, z_step, p, pnew;
double y[3];
x_step = rnd.Gauss(0, a/3);
y_step = rnd.Gauss(0, a/3);
z_step = rnd.Gauss(0, a/3);

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

*/
bool RW :: Accept(double x, double y){
if((x > y)||(rnd.Rannyu() < x/y)) return 1;
else return 0;
}


int RW :: GetCount(void){
return count;
}


double RW :: GetPotential(void){

return x*x*x*x - 2.5*x*x;
}


double RW :: GetIntegrand(void){

return (-0.5*wf.Lap(x) + GetPotential()*wf.Eval(x))/wf.Eval(x);
}



void RW :: Set_mu(double v){
 wf.mu = v;
 return;
}

void RW :: Set_sigma(double v){
 wf.sigma = v;
 return;
}
/*************************************************************************************/
/*************************************************************************************/

/*     Metodi classe Orbital   */

/*************************************************************************************/
/*************************************************************************************/

Wave_func :: Wave_func (){

//mu = 1;
//sigma = 1;
}

Wave_func :: ~Wave_func(){}


double Wave_func :: Lap(double x){                      
 
return Myexp(x,-1)*((x-mu)*(x-mu)/pow(sigma, 4) - 1.0/(sigma*sigma)) + Myexp(x,1)*((x+mu)*(x+mu)/pow(sigma, 4) - 1.0/(sigma*sigma));
}


double Wave_func :: Eval(double x){

 return Myexp(x, -1) + Myexp(x, 1);
}


double Wave_func :: Myexp(double x, int sign){

 return exp(-(x+sign*mu)*(x+sign*mu)/(2.0*sigma*sigma));
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
