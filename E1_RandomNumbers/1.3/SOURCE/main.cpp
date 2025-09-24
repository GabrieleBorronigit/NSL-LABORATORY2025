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
#include <iomanip>
#include <cmath>
using namespace std;




double error(double av, double av2, int n){ //ritorna la deviazione standard della media 

    if(n==0)
    return 0; 

        return sqrt((av2-pow(av,2))/n);

}


 
int main (int argc, char *argv[]){
//initialize the generator
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


//simulating the Buffon experiment, set d=5cm and L=3cm

double d=2.0, L=1.0; //distance between lines and needle length

//naming d= distance between the needle's midpoint (L/2) and the closest line and theta the angle between the needles and the lines,
//one can say that the needle has hit the line when x<=L/2 sin(theta). in the following code, i will sample x as a uniform number in [0,d/2] and 
//a random angle theta using x,y coordinates. Then i will produce an estimation of pi using data blocking method


//pi=2*L*N_hit/d*N_throws
   double x{},y{},sin_theta{};
    int N_throws{1000000};
    int N_blocks{100};
    int Block_size=N_throws/N_blocks;
    string filename="../OUTPUT/Buffon_experiment.dat";
    vector<double> Block_ave;
    vector<double> Block_ave2;
    vector<double> ave;
    vector<double> ave2;
    vector<double> err;
    cerr<<" Throwing needles..."<<endl;
    //loop over blocks
    for(int i{}; i<N_blocks; i++){
      int N_hit{};

    //loop in each block
      for(int j{};j<Block_size;j++){
        //sample a random midpoint position
        double r= rnd.Rannyu(0,d/2);

        //sample sintheta= abs(y/(sqrt(x^2+y^2))) where x,y are r.g. coordinates on the unitary cirle
        do{
         x=rnd.Rannyu(-1,1);
         y=rnd.Rannyu(-1,1);
        }while( pow(x,2)+pow(y,2)==0 || sqrt(pow(x,2)+pow(y,2))>1 ); //to avoid null denominators and points outside the unitary circle
      
       sin_theta=fabs(y/(sqrt(pow(x,2)+pow(y,2))));

       //check hit condition
       if( r<= L/2*sin_theta)
       N_hit++;

      }

      //estimate pi over current block

      double pi= (N_hit>0) ? (2.0*L*Block_size)/(d*N_hit) : 0.0; // if no hits occured, set pi=0
      //update block accumulators
      Block_ave.push_back(pi);
      Block_ave2.push_back(pow(pi,2));
      //loading bar
      if(i%10==0)
      cerr<<i+10<<"% completed"<<endl;
    }

    //loop over blocks
    for(int i{};i<N_blocks; i++){

      double sum_prog{}, sum_prog2{};

      for(int j{}; j<=i; j++){

        sum_prog+=Block_ave[j];
        sum_prog2+=Block_ave2[j];

      }
  
     ave.push_back(sum_prog/(i+1));
     ave2.push_back(sum_prog2/(i+1));
     err.push_back(error(ave[i], ave2[i], i));

    }

    //print results
    ofstream out;
    out.open(filename);

    out<<"#BLOCK"<<setw(12)<<"AVERAGE"<<setw(12)<<"ERROR"<<endl;

    for(int i{};i<N_blocks; i++){

      out<<i+1<<setw(12)<<ave[i]<<setw(12)<<err[i]<<endl;
    }

  out.close();

    


   
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
