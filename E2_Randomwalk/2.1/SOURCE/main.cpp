#include "functions.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "random.h"
#include <string>
using namespace std;

int main() {

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


  // evaluate integral of M_PI/2.0*cos(x*M_PI/2.0) in (0,1) using MC methods
  // uniform distribution
  double a{}, b{1.0};
  int N_throws{10000}, N_blocks{100};
  Integrand f;
  Integral_MC I(a, b, N_throws,rnd);
  vector<double> results;
  cerr<<" evaluating the integral using uniform-sampled random numbers..."<<endl;
  for (int i = 0; i < N_blocks; i++) {
    results.push_back(I.Uniform_MC(f));
    //loading bar
    if(i%(N_throws/10)==0)
    cerr<<(100 * i / N_throws)+10<<"% complete"<<endl;
  }
  // data blocking
  string filename = "../OUTPUT/UniformMC.dat";
  ComputeMeanErrors_DB(results, N_blocks, filename);

  // importance sampling
  // recycle previous vector
  results.clear();
  vector<double> results_linear;

  // allocate PDFs
  IS_P_integrand p;
  IS_P_linear l;      //to show the theoretical "perfect sampling"

  string fn_P_integrand = "../OUTPUT/MCBestIS.dat";
  string fn_P_linear = "../OUTPUT/MCLinearIS.dat";
 cerr<<" evaluating the integral using importance sampling..."<<endl;
  for (int i{}; i < N_blocks; i++) {
    results.push_back(I.Importance_Sampling(f, p));
    results_linear.push_back(I.Importance_Sampling(f, l));
    //loading bar
    if(i%(N_throws/10)==0)
    cerr<<(100 * i / N_throws)+10<<"% complete"<<endl;
  }

  ComputeMeanErrors_DB(results, N_blocks, fn_P_integrand);
  ComputeMeanErrors_DB(results_linear, N_blocks, fn_P_linear);

  cout<<"Done! results are in the files: "<< filename<< "\t"<< fn_P_integrand << "\t"<< fn_P_linear<<endl;



  return 0;
}