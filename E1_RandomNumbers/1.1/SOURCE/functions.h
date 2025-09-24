#pragma once
#include "random.h"
#include <vector>
#include <string>
using namespace std;
double error(double av, double av2, int n);  //deviazione standard della media



// funzione che calcola media e deviazione standard usando data blocking e stampa i risultati su file esterno
void ComputeMeanErrors_DB(const vector<double> &data, int N_blocks, string filename);

//funzione che calcola il chi2
double compute_chi2(int n, const int M, Random &rnd);