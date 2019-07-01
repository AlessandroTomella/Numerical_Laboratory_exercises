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
  Random(int);
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

TSP(int, bool);
TSP(){};
~TSP();
 bool geometry;
vector<double> x;
vector<double> y;

void GenerateCityMap(Random);
void LoadCityMap();
void PrintCityMap();
void SaveCityMap();

int GetCityNum();
virtual double CostFunc(vector<int>);

};


//Percorso del viaggiatore


class Individual: public TSP {

private:
  int size;
  vector<int> ID;
  double fitness;

public:
  //Costruttore	
  Individual(){};	
  Individual(TSP, int, bool);
  //Distruttore
  ~Individual();
  //Accesso alle variabili private
 
 
  void SetElement(int, int);
  int GetElement(int);
  void SetNewID(vector<int>);
  double GetFitness(){return fitness;};
  void ShowInd();
  void PrintInd(int);
  bool SelfCheck();
  bool Check(vector<int>);
  double CostFunc();
};


//Simulated annealing

class SimAnneal : public Random{

private:

  double p_mute;
  int n_shuffle;
  double acc;

public:
  //constructor
  SimAnneal(int, double, int);				//Nel costruttore inzializzo le variabili, tra cui la posizione iniziale
  //destructor
  ~SimAnneal();
  
  double GetAccepted(){return acc;};
  void ResetAccepted(){acc = 0;};
  
  void PairPermutation(Individual &);
  void Shuffle(Individual &);
  void Shift(Individual &);
  void GroupPermutation(Individual &);
  void Inversion(Individual &);
  void RandomSearch(Individual &);
  
  void Metropolis(Individual &, double);
  double Boltzmann(Individual , double);
  //vector<Individual>  MAIN(vector<Individual>); 
  
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
