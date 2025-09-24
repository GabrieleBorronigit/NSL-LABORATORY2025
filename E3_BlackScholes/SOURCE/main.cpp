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

double max(double a, double b){
   return (a>b)? a: b;
}
double error(double av, double av2,
             int n) {  // ritorna la deviazione standard della media

  if (n == 0)
    return 0;  // controllo che non passi dataset vuoti, per evitare divisioni
               // per zero

  return sqrt((av2 - pow(av, 2)) / n);
}


//computes statistical averages and errors using data blocking, outputs on file
void ComputeMeanErrors_DB(const vector<double> &data, int N_blocks, string filename) {

  int N_throws = data.size();

  int Blocks_size = static_cast<int>(N_throws / N_blocks);

  double temp_sum{};
  vector<double> average;
  vector<double> square_average;

  //loop over blocks
  for (int i{}; i < N_blocks; i++) {
    temp_sum = 0;
  //loop inside each block
      for (int j{}; j < Blocks_size; j++) {
        temp_sum += data[i * Blocks_size + j];
       }

    average.push_back(temp_sum / Blocks_size);    // average
    square_average.push_back(pow(average[i], 2)); //square average
  }

 //accumulators
  vector<double> mean_prog(N_blocks,0.0);
  vector<double> sqmean_prog(N_blocks, 0.0);
  vector<double> err_prog(N_blocks,0.0);

  for (int i{}; i < N_blocks; i++) {

    for (int j{}; j <= i; j++) { 
      mean_prog[i] += average[j];
      sqmean_prog[i] += square_average[j];
    }
    mean_prog[i]/=(i+1);      
    sqmean_prog[i]/= (i + 1); 
  
      err_prog[i]=error(mean_prog[i], sqmean_prog[i], i); //<x^2> - <x>^2
  }


  // open file
  ofstream out;
  out.open(filename);
  // check 
  if (!out) {
    cerr << "impossibile aprire il file " << filename << endl;
    return;
  }
  //header
  out<<"#BLOCK \t \t AVERAGE \t \t ERROR"<<endl;
  // print results
  for (int i{}; i < N_blocks; i++) {
    out <<i+1<< " \t \t "<< mean_prog[i] << " \t \t " << err_prog[i] << endl;
  }
  out.close();}




int main (int argc, char *argv[]){

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

   rnd.SaveSeed();

//3.1


double S_0{100.0}, T{1.0}, k{100.0}, r{0.1}, sigma{0.25}, zi{}, S_t{};
int N{100000}, N_blocks{100};
vector<double> C_i, P_i;

//discreto
 for(int i{};i<N; i++)
{
    zi=rnd.Gauss(0,1);
   S_t=S_0*exp((r-pow(sigma,2)/2)*T+sigma*zi*sqrt(T));
   C_i.push_back(exp(-r*T)*max(0,S_t-k));
P_i.push_back(exp(-r*T) * max(k - S_t, 0.0));

}


string filename="../OUTPUT/Call.dat";
string filename2="../OUTPUT/Put.dat";
ComputeMeanErrors_DB(C_i, N_blocks, filename);
ComputeMeanErrors_DB(P_i, N_blocks, filename2);
C_i.clear();
P_i.clear();
double S_i{};
int M=static_cast<int> (100/T);
double dt=T/M;
for(int i{};i<N;i++){
   S_i=S_0;
   for(int j{};j<M;j++)
   {
       zi=rnd.Gauss(0,1);
   S_i=S_i*exp((r-pow(sigma,2)/2)*dt+sigma*zi*sqrt(dt));

   }
   S_t=S_i;
   C_i.push_back(exp(-r*T)*max(0,S_t-k));
  P_i.push_back(exp(-r*T) * max(k - S_t, 0.0));


}

string filename3="../OUTPUT/Call_discrete.dat";
ComputeMeanErrors_DB(C_i, N_blocks, filename3);
string filename4="../OUTPUT/Put_discrete.dat";
ComputeMeanErrors_DB(P_i, N_blocks, filename4);
   return 0;
}


