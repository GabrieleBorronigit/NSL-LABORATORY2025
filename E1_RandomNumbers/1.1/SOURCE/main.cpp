#include "functions.h"
#include "random.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
  Random rnd;
  // initialize generator
  int seed[4];
  int p1, p2;
  ifstream Primes("../INPUT/Primes");
  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else
    cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("../INPUT/seed.in");
  string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;

  // es. 1.1.1 valor medio e incertezza sui blocchi

  int N_throws = 100000;
  int N_blocks = 100;
  vector<double> data;

  // riempio il vettore con numeri random
  for (int i{}; i < N_throws; i++) {
    data.push_back(rnd.Rannyu());
  }

  // calcolo media ed errore sui blocchi
  string fn_testAV = "../OUTPUT/RGtestAV.dat";
  ComputeMeanErrors_DB(data, N_blocks, fn_testAV);

  // 1.1.2 varianza

  // azzero il vettore di prima per riutilizzarlo
  data.clear();

  // riempio il vettore
  for (int i{}; i < N_throws; i++) {
    data.push_back(pow((rnd.Rannyu() - 0.5), 2));
  }
  string fn_testSigma2 = "../OUTPUT/RGtestSIGMA2.dat";
  ComputeMeanErrors_DB(data, N_blocks, fn_testSigma2);

  // 1.1.3
  int M = 100;   // numero di sottointervalli in cui dividere l'intervallo [0,1)
  int n = 10000; // numero di lanci
  double chi[100]; // array che salva i risultati

  // calcolo il chi2 per ogni blocco
  for (int i = 0; i < 100; i++) {

    chi[i] = compute_chi2(n, M, rnd);
  }

  string filename3 = "../OUTPUT/Chi2.dat";
  ofstream out_chi;
  out_chi.open(filename3);

  for (int i{}; i < 100; i++) {
    out_chi << chi[i] << endl;
  }
  out_chi.close();

  return 0;
}