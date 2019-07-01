/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
const double PI = 3.14159265358979323846;
const int EQ_STEPS = 500;			//Numero step nella fase di equilibrazione

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


//Funzione d'onda e derivate

class Wave_func {

private:

public:
  //Variabili (dichiarate pubbliche per accessibilità della classe derivata RW) 
  double mu, sigma;
  
  //Costruttore		
  Wave_func();
  //Distruttore
  ~Wave_func();
 
  double Lap(double);	//Polinomio di Laguerre
  double Myexp(double, int);	//Armonica sferica
  double Eval(double); 	//Modulo quadro funzione d'onda
};


//Random walk 1D con integrato algoritmo di Metropolis

class RW : public Wave_func, public Random{

private:
  Wave_func wf;			//Variabile di tipo Orbitale (nel costruttore introduco i numeri quantici)
  Random rnd;			//Varibile di tipo Random
  double x;			//Posizione nello spazio
  double a;			//Passo del Random Walk
  int count;			//Contatore mosse accettate
public:
  //constructor
  RW();				//Nel costruttore inzializzo le variabili, tra cui la posizione iniziale
  //destructor
  ~RW();
  void Set_mu(double);
  void Set_sigma(double);  
  void ResetCount(void);	
  void SetPoint(double);  
  void Equilibrate(void);	//Fase di equilibrazione
  void PrintPoint(void);

  void Mrt_Step(void);		//Propone la mossa e la accetta secondo l'algoritmo di Metropolis con prob di scelta uniforme
  //void G_Step(void);		//Propone la mossa e la accetta secondo l'algoritmodi Metropolis con prob di scelta Gaussiana
  bool Accept(double, double);	//Accettazione della mossa
  double GetPoint(void);		//Stampa su file la posizione attuale
  int GetCount(void); 		//Mostra variabile count
  double GetPotential(void);	//Calcola il valore dell'energia potenziale
  double GetIntegrand(void);

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
