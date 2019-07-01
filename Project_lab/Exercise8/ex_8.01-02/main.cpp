#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include "funz.h"

//constants	//Numero di cammini generati intotale

const int nblocks = 100;		//Numero di blocchi 
const int dimblocks = 10000;		//Numero step in un blocco

const int pstep = 500;			//Numero tentativi cambio parametri
const double delta = 0.02;

using namespace std;

int main(){  
/*
 int filenum;
 ifstream input;
 input.open("txt");
 input >> filenum;
 input.close();

 cout << "ESECUZIONE NUMERO " << filenum << endl << endl;
*/ 
 int conf = 1;
 double E, ave_E = 0, var_E = 0;
 double mean = 0.802939, sigma = 0.620538;
/* if (filenum < 10){
 mean = (filenum) * 1.0/10.0;
 sigma = 1;
 }
 else {
 sigma = (filenum - 9) * 1.0/10.0;
 mean = 1;
 }*/
 RW walker;
 walker.Set_mu(mean);
 walker.Set_sigma(sigma);

walker.Equilibrate();
 cout << "Starting the sampling at the ";
 walker.PrintPoint();
 walker.ResetCount(); 
 
/* 
 for(int i = 1; i <= nblocks; i++){		
 
	E = 0;

	for(int j = 1; j <= dimblocks; j++){
 	
	E += walker.GetIntegrand()/(double)dimblocks;
 	walker.Mrt_Step();	
	
	}
	
	ave_E += E;

 }

 cout << "Starting value for the energy (MU = "<< mean << ", SIGMA = "<< sigma << "): " << ave_E/(double)nblocks << endl;
 
 ofstream param, energy;
 param.open("HDwalk.dat");
 //param.open("parameters_walk/walk."+ to_string(filenum+56) + ".dat");
 param << mean << " " << sigma << endl;
 
 int mod = 0;
 double trial_mean, trial_sigma, ave_Ebis = 0;
 walker.ResetCount();
 double step = delta;
 
 for(int k = 1; k <= pstep; k++){

 if(k > 100) step = delta/2.0;
 if(k%2 == 0){ trial_mean = mean + step*(walker.Rannyu()-0.5);
 	walker.Set_mu(trial_mean);
 }
 if(k%2 == 1){trial_sigma = sigma + step*(walker.Rannyu()-0.5); 
 walker.Set_sigma(trial_sigma);
 }

 for(int i = 1; i <= nblocks; i++){		
 
	E = 0;

	for(int j = 1; j <= dimblocks; j++){
 	
	E += walker.GetIntegrand()/(double)dimblocks;
	
 	walker.Mrt_Step();	
	
	}
	
	ave_Ebis += E;
  }

  cout << k << ") Accepted to proposed move ratio: " << (double)walker.GetCount()/(double)(nblocks*dimblocks) << endl;
  
  if(ave_Ebis < ave_E){
 	mod++;
	mean = trial_mean;
  	sigma = trial_sigma;
  	ave_E = ave_Ebis;
	
   }

 param << mean << " " << sigma << endl;
 
 ave_Ebis = 0;
 
 walker.ResetCount(); 

}
 
 cout << "Final Parameters: " << mean << " " << sigma << endl;
 cout << "Lowest energy found: " << ave_E/(double)nblocks << endl;
 cout << "Accepted parameter changes:" << mod << endl;
 param.close();	

 if(ave_E/(double)nblocks > -0.4){
 cout <<"//////////////BAD MINIMIZATION/////////////"<< endl << endl;
 return 0; 
 }
 else{

	energy.open("Energy.dat", ios::app);
	energy << filenum+56 << " " << ave_E/(double)nblocks << endl;
 	energy.close();

 }
 
*/

 //cout << "Saving optimal-value parameters mu = " << mean <<", sigma = " << sigma << endl << endl;
 //walker.Set_mu(mean);
 //walker.Set_sigma(sigma);


 ofstream points, Energy;
 points.open("sample.out");
 Energy.open("Ave_en.out");
 
 ave_E = 0;

 for(int i = 1; i <= nblocks; i++){		
 
	E = 0;
	cout << "Blocco numero " << i << endl;
	for(int j = 1; j <= dimblocks; j++){
 	
	E += walker.GetIntegrand()/(double)dimblocks;
 	walker.Mrt_Step();	
	points << walker.GetPoint() << endl;
	}
	
	ave_E += E;
	var_E += E*E;
	
	Energy << fixed << setprecision(7) << ave_E/(double)i << " " << sqrt((var_E/(double)i - pow(ave_E/(double)i, 2))/(double)i)<< endl;

 }

 points.close();
 Energy.close();


return 0;
}

