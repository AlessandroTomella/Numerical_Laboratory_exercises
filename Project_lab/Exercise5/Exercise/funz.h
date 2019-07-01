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


//Orbitale atomico

class Orbital {

private:
  //Numeri quantici
  int n, m, l;

public:
  //Costruttore		
  Orbital();
  //Distruttore
  ~Orbital();
  //Accesso alle variabili private
  int get_n();
  int get_l();
  int get_m();
  
  int Fact(int);		//Fattoriale
  double Bin(int, int);		//Coefficiente Binomiale
  double Norm(void);		//Normalizzazione orbitale
  double Laguerre(double);	//Polinomio di Laguerre
  double Harmonic(double);	//Armonica sferica
  double Eval(double r[3]); 	//Modulo quadro funzione d'onda
};


//Random walk con integrato algoritmo di Metropolis

class RW{

private:
  Orbital orb;			//Variabile di tipo Orbitale (nel costruttore introduco i numeri quantici)
  Random rnd;			//Varibile di tipo Random
  double x[3];			//Posizione nello spazio
  double a;			//Passo del Random Walk
  int count;			//Contatore mosse accettate
  int EQ_STEPS;			//Numero step nella fase di equilibrazione
public:
  //constructor
  RW();				//Nel costruttore inzializzo tutte le variabili, tra cui la posizione iniziale
  //destructor
  ~RW();
  
  void Equilibrate(void);	//Fase di equilibrazione
  void PrintPoint(void);	//Stampa a video la posizione corrente
  void U_Step(void);		//Propone la mossa e la accetta secondo l'algoritmo di Metropolis con prob di scelta uniforme
  void G_Step(void);		//Propone la mossa e la accetta secondo l'algoritmodi Metropolis con prob di scelta Gaussiana
  bool Accept(double, double);	//Accettazione della mossa
  double GetCoord(int);		//Stampa su file le coordinate del punto attuale
  int GetCount(void); 		//Mostra variabile count
  double GetRadius(void);	//Calcola la distanza dall'origine della posizione attuale

  int Getn(void);		//Funzioni di accesso ai numeri quantici della variabile orb
  int Getl(void);
  int Getm(void);
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
