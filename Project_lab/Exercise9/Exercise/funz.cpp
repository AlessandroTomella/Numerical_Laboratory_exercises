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

Random :: ~Random(){cout << "***" << endl;}//SaveSeed();}

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


void TSP :: GenerateCityMap(Random gen){
 // Generate cities on circle
 if(geometry){
	vector<double> theta;

	for (int i = 0; i < citynum; i++) theta.push_back(gen.Rannyu(0, 2*PI));
	sort(theta.begin(), theta.end());
	for (int i = 0; i < citynum; i++) {
		x.push_back(cos(theta.at(i)));
		y.push_back(sin(theta.at(i)));
	}
 }
 // Generate cities inside a square
 else{
	for (int i = 0; i < citynum; i++) {
		x.push_back( gen.Rannyu());
		y.push_back( gen.Rannyu());
	}
 }

return;
}

// Load cities from file
void TSP :: LoadCityMap(){
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

// Return cost of path v
double TSP :: CostFunc(vector<int> v){

double L = 0;
for(int i = 0; i < v.size()-1; i++) {

	L += pow(x.at(v.at(i)) - x.at(v.at(i+1)), 2) + pow(y.at(v.at(i)) - y.at(v.at(i+1)),2);
}

return L + pow(x.at(v.at(0)) - x.at(v.at(v.size()-1)), 2) + pow(y.at(v.at(0)) - y.at(v.at(v.size()-1)),2);

}

/*************************************************************************************/
/*************************************************************************************/

/*     Metodi classe GenAlg   */

/*************************************************************************************/
/*************************************************************************************/


GenAlg :: GenAlg(int init_shuffle, double pselection, double pmutation, double pcrossover){

n_shuffle = init_shuffle;
roulette = pselection;
p_mute = pmutation;
p_cross = pcrossover;


}

GenAlg :: ~GenAlg(){SaveSeed();}

// This function sort the Population from low to high cost
// This could be done also by overloading the "<" operator
void GenAlg :: Hierarchy(vector<Individual> & Population){

 if(Population.size() == 1) {
 cout << "ERROR: Cannot order a size=1 Array" << endl;
 return;
 }
 vector<double> OriginalFitness, OrderedFitness;
 for(int i = 0; i< Population.size(); i++) {
 	OriginalFitness.push_back(Population.at(i).CostFunc());
	OrderedFitness.push_back(OriginalFitness.at(i));
 }
 // Sort a vector of doubles
 sort(OrderedFitness.begin(), OrderedFitness.end());

 vector<Individual> CopyPopulation;
 // Reconstruct the sorted population from the sorted vector
 for(int i = 0; i < OriginalFitness.size(); i++){
 	for(int k = 0; k < OriginalFitness.size(); k++){
		if(OriginalFitness.at(k) == OrderedFitness.at(i)){
			CopyPopulation.push_back(Population.at(k));
			break;
		}		
 	}
 }
 // Set the old disordered population equal to the sorted one
 Population = CopyPopulation;

return;
}


// Permuta due elementi del DNA dell'individuo
void GenAlg :: PairPermutation(Individual & cavy){

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

// Esegue una permutazione multipla sul DNA dell'individuo
void GenAlg :: Shuffle(Individual & cavy){

 for(int i = 0; i < n_shuffle; i++){
 PairPermutation(cavy);
 }

return;
}


// Shifta una zona del DNA dell'individuo in un'altra zona 
void GenAlg :: Shift(Individual & cavy){

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


// Scambia due pezzi del DNA dell'individuo
void GenAlg :: GroupPermutation(Individual & cavy){

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


// Esegue un'inversione di una zona di DNA dell'individuo
void GenAlg :: Inversion(Individual & cavy){

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

// La randomsearch chiama con una certa probabilità una delle alterazioni genetiche appena elencate
void GenAlg :: RandomSearch(Individual & cavy){

 double r = Rannyu();

 if(r < p_mute) PairPermutation(cavy);
 	
 if(r > p_mute && r < 2*p_mute) Shift(cavy);
	
 if(r > 2*p_mute && r < 3*p_mute) GroupPermutation(cavy);
	
 if(r > 3*p_mute && r < 4*p_mute) Inversion(cavy);
 
 return;
}


// Scambia elementi di DNA tra due individui
void GenAlg :: CrossOver(Individual & Father, Individual & Mother){

 vector<int> Son1, Son2;

 int dim = Mother.GetCityNum();
 int Cut_ind = (int) Rannyu(1, dim-2);		// punto in cui effettuare il taglio

 for(int i = 0; i < Cut_ind; i++){
	Son1.push_back(Mother.GetElement(i));	// copia nei figli la prima parte del DNA 
	Son2.push_back(Father.GetElement(i));
 }   

 vector<int> FCutted(dim - Cut_ind);
 vector<int> FPosition(dim - Cut_ind);
 vector<int> MCutted(dim - Cut_ind);
 vector<int> MPosition(dim - Cut_ind);

 for(int i= 0; i< dim - Cut_ind; i++){
	FCutted[i] = Father.GetElement(i+Cut_ind);	// Salva i pezzi tagliati in nuovi vettori
	MCutted[i] = Mother.GetElement(i+Cut_ind);
 }

 for (int i = 0; i < dim - Cut_ind; i++){				// "Se l'elemento j-esimo del genitore1 è
	for (int j = 0; j < dim ; j++){					//  uguale all'elemento i-esimo della parte tagliata
		if(Mother.GetElement(j) == FCutted[i]) MPosition[i] = j;//  al genitore2, allora salva la posizione j all'i-esimo posto
		if(Father.GetElement(j) == MCutted[i]) FPosition[i] = j;//  di un nuovo vettore *Position[]
	}
 }

 sort(FPosition.begin(),FPosition.end());	
 sort(MPosition.begin(),MPosition.end());

 
 for(int i = 0; i < dim - Cut_ind; i++){
	Son1.push_back(Father.GetElement(FPosition[i]));
	Son2.push_back(Mother.GetElement(MPosition[i]));
 }   

 Mother.SetNewID(Son1);
 Father.SetNewID(Son2);

return ; 
}


// Gestisce un passo evolutivo e restistuisce la nuova popolazione
vector<Individual>  GenAlg :: Evolve(vector<Individual> Population){
 int dim = Population.size(), count = 0;
 vector<Individual> New_Population;	// Nuova popolazione
 double pos1, pos2;

 while(New_Population.size() < dim){ 	// ripetere finchè la dimensione della nuova popolazione è uguale a quella vecchia

	pos1 = pow(Rannyu(), roulette)*dim;	// Selezione a roulette truccata che favorisce 
	pos2 = pow(Rannyu(), roulette)*dim; 	// i primi elementi della popolazione (già ordinata in questa fase)

	for(int i = 0; i < dim; i++){		// Creazione a 2 a 2 della nuova popolazione
		if(pos1 > i && pos1 < i+1) New_Population.push_back(Population.at(dim - i - 1));
		if(pos2 > i && pos2 < i+1) New_Population.push_back(Population.at(dim - i - 1));
	}

	if(Rannyu()< p_cross) CrossOver(New_Population.at(count), New_Population.at(count+1));	// Crossover con probabilità p_cross

	for (int i = count; i < count + 2; i++){
	RandomSearch(New_Population.at(i));		// Algoritmo di Randomsearch
	}
count+= 2;	// Al termine aggiorna il contatore di individui di 2
}

for(int i = 0; i< New_Population.size(); i++){		// Controlla quello che è stato creato
if(!New_Population.at(i).SelfCheck()) cout << "***ERROR***" << endl;
}

return New_Population;	// Restituisce la nuova popolazione

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
 if(geometry) print.open("WorkInProg/Circular/CirclePath."+to_string(file)+".final");
 else print.open("WorkInProg/Square/SquarePath."+to_string(file)+".final");
 
 for(int i = 0; i < size; i++) {	
     	   print << ID.at(i) <<  endl;
 	}
 
 print.close();

return;
}


bool Individual :: Check(vector<int> Seq){

 double x;
 
 if (Seq.size() != 30) {
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


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
