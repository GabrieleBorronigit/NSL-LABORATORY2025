#include "Genetic_algorithm.h"
#include "random.h"
#include <armadillo>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string>
using namespace std;
using namespace arma;

int main(int argc, char *argv[]) {

  // generate the population
  int type= stoi(argv[1]);
  if(type!=0 && type!=1){
    cerr<<"type 0 for cities randomly placed on a circumference and 1 for cities randomply placed on a square"<<endl;
    exit(1);
  }

  int N_cities = 34;
  int N_chromosomes=200;
  int N_generations = 3500;

  mat cities = Generate_cities(N_cities, type);
  // compute distances
  mat D = Distances(cities);

  // genetic algorithm

  GA Darwin( N_chromosomes, N_cities);

  ofstream out;
  out.open("../OUTPUT/besthalfloss.dat");
  ofstream outbest("../OUTPUT/Best_loss.dat");
  out<<"#gen      BHL"<<endl;
  for (int i = 0; i < N_generations; i++) {

    Darwin.Evolve(D);
    
    Darwin.Sort_by_fitness(D);
     out<<setw(12)<<i<<setw(12)<<Darwin.BestHalfLoss(D)<<endl;
    outbest<<setw(12)<<i<<setw(12)<<Darwin.Loss(D,0)<<endl;
  }
  out.close();
  outbest.close();

  Darwin.Pick_best(D);

  return 0;
}
