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
#include <bits/stdc++.h> 
#include <iomanip>
#include "vector"
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
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            SetRandom(seed,p1,p2);

	input.close();
      }
      
    else cerr << "PROBLEM: Unable to open seed.in" << endl;
return;
}

Random :: ~Random(){}

void Random :: SaveSeed(){
//cout << "Saving seed" << endl;
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
//cout << "Numero random: " <<r << endl;
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

/*     Metodi classe TPS   */

/*************************************************************************************/
/*************************************************************************************/

TSP :: TSP(int ncity, bool geom){
 citynum = ncity;
 geometry = geom;
}

TSP :: ~TSP(){} 

void TSP :: GenerateCityMap(Random ptr){

cout << "Generating a new random distribution of cities" << endl; 

 if(geometry){
	vector<double> theta;

	for (int i = 0; i < citynum; i++) theta.push_back(ptr.Rannyu(0, 2*PI));
	sort(theta.begin(), theta.end());
	for (int i = 0; i < citynum; i++) {
		x.push_back(cos(theta.at(i)));
		y.push_back(sin(theta.at(i)));
	}
 }

 else{
	for (int i = 0; i < citynum; i++) {
		x.push_back( ptr.Rannyu());
		y.push_back( ptr.Rannyu());
	}
 }

return;
}


void TSP :: LoadCityMap(){

cout << "Loading a previous distribution of cities" << endl; 
 double a, b;
 ifstream input;
 if(geometry) input.open("circle.start");
 else input.open("square.start");
 
 if(input){
 	for (int i = 0; i < citynum; i++) {
	input >> a >> b;
	x.push_back(a);
	y.push_back(b);
	}
 }
 else cout << "Map of cities not found " << endl;
 input.close();

return;
}



void TSP :: SaveCityMap(){

cout << "Saving the newly generated distribution of cities" << endl; 

 ofstream print;
 if(geometry) print.open("circle.start");
 else print.open("square.start");
 for (int i = 0; i < citynum; i++) print <<  fixed << setprecision(7) << x[i] << " " << y[i] << endl;
 
 print.close();
}


void TSP :: PrintCityMap(){
 cout << "Coordinate delle città sulla mappa: " << endl;
 cout << "    x    " << "    y    " << endl;
 for (int i = 0; i < citynum; i++) cout <<  fixed << setprecision(7) << x[i] << " " << y[i] << endl;
 
}


int TSP :: GetCityNum(){

return citynum;

}


double TSP :: CostFunc(vector<int> v){

double L = 0;

for(int i = 0; i < v.size()-1; i++) {

	L += pow(x.at(v.at(i)) - x.at(v.at(i+1)), 2) + pow(y.at(v.at(i)) - y.at(v.at(i+1)),2);
}

return L + pow(x.at(v.at(0)) - x.at(v.at(v.size()-1)), 2) + pow(y.at(v.at(0)) - y.at(v.at(v.size()-1)),2);

}

/*************************************************************************************/
/*************************************************************************************/

/*     Metodi classe SimAnneal   */

/*************************************************************************************/
/*************************************************************************************/


SimAnneal :: SimAnneal(int init_shuffle, double pmutation){

n_shuffle = init_shuffle;
p_mute = pmutation;
acc = 0;
}

SimAnneal :: ~SimAnneal(){SaveSeed();}

void SimAnneal :: PairPermutation(Individual & cavy){

 int rand1, rand2, temp;
 int range = cavy.GetCityNum();
 rand1 = (int) Rannyu(0, range);
 rand2 = (int) Rannyu(0, range);
 while(rand1 == rand2) rand2 = (int) Rannyu(0, range);

 temp = cavy.GetElement(rand1);
 cavy.SetElement(cavy.GetElement(rand2), rand1);
 cavy.SetElement(temp, rand2); 

return;
}


void SimAnneal :: Shuffle(Individual & cavy){

 for(int i = 0; i < n_shuffle; i++){
 PairPermutation(cavy);
 }

return;
}

void SimAnneal :: Shift(Individual & cavy){

 int dim = cavy.GetCityNum();
 int K_shift = (int)Rannyu(1, dim);
 vector<int> temp(dim);
 
 for (int i = 0; i < dim;i++){
 temp[i] = cavy.GetElement(i - K_shift);
 }
 
 for (int i = 0; i < dim;i++){
 cavy.SetElement(temp[i], i);
 }
 
return;
}


void SimAnneal :: GroupPermutation(Individual & cavy){

 int dim = cavy.GetCityNum();
 int M_perm = (int) Rannyu(2, (int) dim/2.0);
 int max_ind = (int) dim/(double)M_perm;
 int first_blk = (int) Rannyu(0, max_ind);
 int sec_blk = (int) Rannyu(0, max_ind);
 while(sec_blk == first_blk) sec_blk = (int) Rannyu(0, max_ind);
 
 vector<int> temp1(M_perm);
 vector<int> temp2(M_perm);

 for (int i = 0; i < M_perm; i++){
 temp1[i] = cavy.GetElement(i + first_blk*M_perm);
 temp2[i] = cavy.GetElement(i + sec_blk*M_perm);
 }
 
 for (int i = 0; i < M_perm; i++){
 cavy.SetElement(temp2[i], i+first_blk*M_perm) ;
 cavy.SetElement(temp1[i], i+sec_blk*M_perm) ;
 }
 
return;
}



void SimAnneal :: Inversion(Individual & cavy){

 int dim = cavy.GetCityNum();
 int M_inv = (int) Rannyu(2, dim-1);
 int start_ind = (int) Rannyu(0, dim);
 vector<int> temp(M_inv);

 for (int i = 0; i < M_inv; i++){
 temp[i] = cavy.GetElement(i + start_ind);
 }
 
 for (int i = 0; i < M_inv;i++){
 cavy.SetElement(temp[M_inv - i - 1], i+start_ind);
 }
 
return;
}


void SimAnneal :: RandomSearch(Individual & cavy){

 double r = Rannyu();

 if(r < p_mute){
	//cout << "Performing a permutation of pair" << endl; 
	PairPermutation(cavy);
 	}
 //r = Rannyu();
 if(r > p_mute && r < 2*p_mute){ 
	//cout << "Performing a shift" << endl;	
	Shift(cavy);
	}
 //r = Rannyu();
 if(r > 2*p_mute && r < 3*p_mute){ 
	//cout << "Performing a group permutation" << endl;
	GroupPermutation(cavy);
	}
 //r = Rannyu();
 if(r > 3*p_mute && r < 4*p_mute){ 
	//cout << "Performing an inversion" << endl;
	Inversion(cavy);
	}
 

 return;
}


double SimAnneal :: Boltzmann(Individual ind, double beta){
double p;

p = exp(-beta*ind.CostFunc());
return p;
}

void SimAnneal :: Metropolis(Individual & ind, double beta){

Individual copy = ind;
RandomSearch(copy);
if (!copy.SelfCheck()) return;
if(Rannyu() < Boltzmann(copy, beta)/Boltzmann(ind, beta)) {
	ind = copy;
	acc++;
}
 
return ; 
}



/*************************************************************************************/
/*************************************************************************************/

/*     Metodi classe Individual   */

/*************************************************************************************/
/*************************************************************************************/

Individual :: Individual(TSP Realization, int ncities, bool geometry) : TSP(ncities, geometry){

 size = Realization.GetCityNum();
 for(int i = 0; i < size; i++){
 x.push_back(Realization.x.at(i));
 y.push_back(Realization.y.at(i));
 ID.push_back(i);
 }
}

Individual :: ~Individual(){}

void Individual :: ShowInd(){

//Stampa il vettore
cout << "Sequenza città: " ;
 for(int i = 0; i < size; i++) {
	 cout << ID.at(i) << "-";
 }
 cout << endl;
 return;
}


void Individual :: PrintInd(int file){

 ofstream print;
 if(geometry) print.open("Circle/path."+to_string(file)+".dat");
 else print.open("Square/path."+to_string(file)+".dat");
 
 for(int i = 0; i < size; i++) {	
     	   print << ID.at(i) <<  endl;
 	}
 
 print.close();

return;
}


bool Individual :: Check(vector<int> Seq){

 double x;
 
 if (Seq.size() != size) {
 cout << "Check Function Error : !SIZE! " << endl;
 return 0;
 }

 for (int i = 0; i < size-1; i++){
 x = Seq.at(i);
	for (int j = i+1; j < size; j++){
 	if(x==Seq.at(j)) {
		cout << "Check Function Error : !CORRUPTED INDIVIDUAL WITH EQUAL ENTRIES!" << endl;		
		return 0;
		}
	} 
 }

 return 1;

}

bool Individual :: SelfCheck(){
 
 return Check(ID);
}


void Individual :: SetElement(int newid, int pos){

 while(pos >= size || pos < 0) {
	if(pos >= size) pos = pos - size;
	if(pos < 0) pos = pos + size;
 }
 ID.at(pos) = newid;

 return;
}

int Individual :: GetElement(int pos){
 
 while(pos >= size || pos < 0) {
	if(pos >= size) pos = pos - size;
	if(pos < 0) pos = pos + size;
 }

 return ID.at(pos);
}

void Individual :: SetNewID(vector<int> NewID){
 
 if(!Check(NewID)){
 cout << "ERROR: Trying to set corrupted ID in sane Individual." << endl;
 return;
 }
 
 for(int i = 0; i < size; i++){
 ID.at(i) = NewID.at(i);
 }

 return;
}


double Individual :: CostFunc(){

 fitness = TSP::CostFunc(ID);
 return fitness;

}

double Individual :: Mod_CostFunc(vector<double> x, vector<double> y){

double L = 0;
for(int i = 0; i < size-1; i++) {

	L += pow(x.at(ID.at(i)) - x.at(ID.at(i+1)), 2) + pow(y.at(ID.at(i)) - y.at(ID.at(i+1)),2);
}

return L + pow(x.at(ID.at(0)) - x.at(ID.at(size-1)), 2) + pow(y.at(ID.at(0)) - y.at(ID.at(size-1)),2);

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
