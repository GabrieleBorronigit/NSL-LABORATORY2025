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
#include <string>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){
//inizializzo il generatore
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("../INPUT/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../INPUT/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  //10000 stime di S_N=1/N sum(x_i) per N=(1,2,10,100),
  // con xi generato random lanciando un dado standard, esponenziale e lorentziano

   int samples=10000;
   int Ns[4]={1,2,10,100};                     //iterazioni
   double S_N_Std{}, S_N_Exp{}, S_N_Lorentz{}; //accumulatori per tipo di dado
   string fnStd="../OUTPUT/Standard_dice.dat", fnExp="../OUTPUT/Exp_dice.dat", fnLor="../OUTPUT/Lorentzian_dice.dat";
   ofstream outstd, outexp, outlor;
   outstd.open(fnStd);
   outexp.open(fnExp);
   outlor.open(fnLor);

    for(int N:Ns)
    {
      //riempio gli accumulatori
      for(int j=0;j<samples;j++)
      {
         S_N_Std=0;
         S_N_Exp=0;
         S_N_Lorentz=0;
         for(int k=0; k<N;k++)
         {
             S_N_Std+=rnd.Rannyu();
             S_N_Exp+=rnd.Exp();
             S_N_Lorentz+=rnd.Lorentz();
         }
         outstd<<N<<setw(16)<<S_N_Std/N<<endl;
         outexp<<N<<setw(16)<<S_N_Exp/N<<endl;
         outlor<<N<<setw(16)<<S_N_Lorentz/N<<endl;
      }

    }
    
outstd.close();
outexp.close();
outlor.close();
cout<<"i dati sono stati esportati nei file \t"<<fnStd<<"\t"<<fnExp<<"\t"<<fnLor<<endl;
   rnd.SaveSeed();
   return 0;
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
