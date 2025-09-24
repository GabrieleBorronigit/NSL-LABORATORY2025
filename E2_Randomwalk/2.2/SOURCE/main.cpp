#include "RandomWalk.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
using namespace std;

int main(){

    Random rnd;
   //initialize generator
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




//set RW parameters
    int N_steps{100};
    int N_walkers{10000};
    int N_blocks{100};
    double a{1.0};



//save two RWs to plot images
RandomWalk RW(N_steps, a, rnd);
bool save_history= true;
RW.LatticeWalk(save_history);
string fn="../OUTPUT/completeRW_lattice.dat";
RW.History_ToFile(fn);

string fn2="../OUTPUT/completeRW_continuum.dat";
RW.Reset();
RW.ContinuumWalk(save_history);
RW.History_ToFile(fn2);
RW.Reset();

//print and accumulate data
vector<vector<double>> amsd_lattice(N_steps, vector<double>(N_walkers, 0.0));
vector<vector<double>> amsd_continuum(N_steps, vector<double>(N_walkers, 0.0)); //matrix that saves, for each step, the distance from origin of every walker
ofstream out_lattice;
ofstream out_continuum;
string fn_lattice="../OUTPUT/LatticeAMSD.dat", fn_continuum="../OUTPUT/ContinuumAMSD.dat";
//check
out_lattice.open(fn_lattice);
out_continuum.open(fn_continuum);
if (!out_lattice) {
    cerr << "Unable to open file " << fn_lattice << endl;
    return 1;
  }
  if (!out_continuum) {
    cerr << "Unable to open file " << fn_continuum << endl;
    return 1;
  }
//headers
 out_lattice<<"#step \t \t AVERAGE \t \t ERROR"<<endl;
 out_continuum<<"#step \t \t AVERAGE \t \t ERROR"<<endl;

 cout<<" walking..."<<endl;

  for(int i{};i<N_walkers;i++){ //setup the walkers

   RW.Reset();
   RW.LatticeWalk(); //walk, for each step, R^2 is saved in a data member of RW
   for(int k{};k<N_steps;k++){
   amsd_lattice[k][i]=RW.Get_amsd_lattice_comp(k); //save every walker's position
   }
   //repeat for continuum
   RW.Reset();
   RW.ContinuumWalk();
    for(int k{};k<N_steps;k++){
   amsd_continuum[k][i]=RW.Get_amsd_lattice_comp(k); 
   
   }
   if(i==N_walkers/2)
   cout<<"collecting AMSD for "<<N_walkers<<" walkers..."<<endl;

  }

//data blocking on all the walkers for every step
for(int i{};i<N_steps;i++){
   ComputeMeanErrors_DB(amsd_lattice[i], N_blocks, out_lattice, i); //prints only the final values of the progressive average
   ComputeMeanErrors_DB(amsd_continuum[i], N_blocks, out_continuum, i);
}

out_lattice.close();
out_continuum.close();

cout<<"Simulation completed! results stored in files:\t"<<fn_lattice<<"\t"<<fn_continuum<<endl;



    return 0;
}