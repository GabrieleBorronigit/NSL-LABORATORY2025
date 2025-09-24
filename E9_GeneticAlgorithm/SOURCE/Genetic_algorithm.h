#pragma once
#include "random.h"
#include <armadillo>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>
using namespace std;
using namespace arma;

class Chromosome {

public:
  Chromosome(int N_genes,
             Random &_rnd); // constructor //passo l indirizzo di memoria del
                            // generatore di GA, in modo da usarne uno per tutti
  Chromosome(); // copy constructor

  bool Check_Boundary(); // boundary condition checker
  void
  Swap(int Ntimes); // swaps two random genes in the chromosome Ntimes times
  const vec &Get_Genes() const { return _genes; }
  double Get_Gene(int i) const { return _genes(i); }
  int Get_Ngenes() const { return _N_genes; }
  void Set_Gene(int i, double gene) { _genes(i) = gene; }
  void Set_Genes(const vec &genes); // sets a whole chromosome (to generate son)
  void
  Swap_genes(int i,
             int j); // swaps the position of two given genes in the chromosome

protected:
  vec _genes;
  int _N_genes;

  Random _rnd;
};

class GA {

public:
  GA(int N_chromosomes, int N_genes);

  // population acces methods
  const Chromosome &Get_Chromosome(int i) { return _population[i]; };
  void Swap_in_population(int i,
                          int j); // swaps two chromosomes in the population
  int GetN();// genetic algorithm

  double Loss(const mat &D, int i); // loss function-to be minimized
  void Sort_by_fitness(const mat &D); // sorts by fitness, starting from the best one
  int Selection(const mat &D);   // selection operator

  // mutation operators
  void Random_mutation(Chromosome &c);
  void Pair_swap(Chromosome &c); // swaps a random pair of cities in the i-chromosome
  void Shift_contiguous(Chromosome &c); // shifts by a random increment a random number of contiguos cities in the i-chromosome
  void Permutation_contiguous(Chromosome &c); // permutation of a random number of contiguous cities in the i-chromosome
  void Invert_path(Chromosome &c); // swaps the positions in the i-chromosomes of a
                              // random number of cities evolution methods
  void Crossover(Chromosome &f1, Chromosome &f2);
  void Evolve(const mat &D); // creates a new population
  void Pick_best(const mat &D); // picks the best chromosome and outputs its genome
  double BestHalfLoss(const mat &D);
  bool TestCrossover(int trials);

protected:
  std::vector<Chromosome> _population;
  Random _rnd;
  int _N_chromosomes;
  double _P_pair, _P_shift, _P_perm, _P_invert; // mutation probabilities
  double _P_crossover;                          // crossover probability
};

// initialization of TSP
mat Generate_cities(int N, int type);
mat Distances(const mat &cities);