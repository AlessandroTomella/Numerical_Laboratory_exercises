/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <vector>
using namespace std;
const double PI = 3.14159265358979323846;

//Generatore di numeri casuali

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  //Random(int *);
  // destructor
  ~Random();
  // methods

  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
};


//Mappa delle città

class TSP {

private:

 int citynum;

public:
//costruttore
TSP(int, bool);
// distruttore
~TSP();
 // se geometry = 1 -> città sul cerchio, altrimenti città nel quadrato
 bool geometry;
// Coordinate delle città
vector<double> x;
vector<double> y;

void GenerateCityMap(Random);
void LoadCityMap();
void PrintCityMap();
void SaveCityMap();

int GetCityNum();
// Funzione virtual perchè l'implementazione 
virtual double CostFunc(vector<int>);

};


//Percorso del viaggiatore


class Individual: public TSP {

private:
// dimensione del cammino
  int size;
// vettore di interi che identifica univocamente il cammino
  vector<int> ID;
// Costo del cammino
  double fitness;

public:
  //Costruttore		
  Individual(TSP, int, bool);
  //Distruttore
  ~Individual();
 
  void SetElement(int, int); 			// cambia elemento dell'ID
  int GetElement(int);				// restituisce elemente dell'ID
  void SetNewID(vector<int>);			// cambia ID
  double GetFitness(){return fitness;};
  void ShowInd();				// stampa l'ID a video
  void PrintInd(int);				// salva l'ID su file
  bool SelfCheck();				// controlla la correttezza dell'ID
  bool Check(vector<int>);
  double CostFunc();				// calcola il costo del cammino
  
};


//Algoritmo genetico

class GenAlg : public Random{

private:
// probabilità di crossover
  double p_cross;
// probabilità mutazioni
  double p_mute;
// esponente della roulette truccata per la selezione
  double roulette;  
// # mescolamenti individui inziali
  int n_shuffle;

public:
  //constructor
  GenAlg(int, double, double, double);
  //destructor
  ~GenAlg();

// Spiegazione delle funzioni nella loro definizione
  void Hierarchy(vector<Individual> &);
  void PairPermutation(Individual &);
  void Shuffle(Individual &);
  void Shift(Individual &);
  void GroupPermutation(Individual &);
  void Inversion(Individual &);
  
  void RandomSearch(Individual &);
  void CrossOver(Individual &, Individual &);

  vector<Individual>  Evolve(vector<Individual>); 
  
};


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
