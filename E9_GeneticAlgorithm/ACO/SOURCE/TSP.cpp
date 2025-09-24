#include "ACO.h"
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
  int type = stoi(argv[1]);
  if (type != 0 && type != 1) {
    cerr << "type 0 for cities randomly placed on a circumference and 1 for "
            "cities randomply placed on a square"
         << endl;
    exit(1);
  }

  int Ncities = 34;
  int Nants = 300;
  int Nwalks = 10;
  double alpha =1.0;
  double beta = 2.0;
  double rho=0.3;
  double Q=100;
  mat cities = Generate_cities(Ncities, type);
  mat D = Distances(cities);
  ofstream outbl("../OUTPUT/best_loss.dat");
  ofstream outbhl("../OUTPUT/besthalfloss.dat");
  Colony anthill(Ncities, Nants, rho, Q);

  for (int i = 0; i < Nwalks; i++) {
    anthill.Search(D, alpha, beta);
    outbl << anthill.Best_loss(D) << endl;
    outbhl << anthill.BestHalfLoss(D) << endl;
  }
  anthill.Pick_best(D);
outbl.close();
outbhl.close();
  return 0;
}