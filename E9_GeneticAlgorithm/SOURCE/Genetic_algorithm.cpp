
#include "Genetic_algorithm.h"
#include "random.h"
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;
using namespace arma;

////////////////////////INITIALIZATION OF TSP METHODS///////////////////////////

mat Generate_cities(int N, int type) {

  Random _rnd;
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
        _rnd.SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;
  _rnd.SaveSeed();

  ofstream outcities;
  outcities.open("../OUTPUT/Cities_Map.dat");
  outcities << "#city               x               y" << endl;
  mat cities(N, 2);

  cities.zeros(N, 2);

  if (type == 0) // randomly placed on a circumference
  {
    // sample sintheta and costheta ( R=10)
    double R = 1.0;
    double theta = 2 * M_PI * _rnd.Rannyu();

    for (int i = 0; i < N; i++) {

      double theta = 2 * M_PI * _rnd.Rannyu();
      cities(i, 0) = R * cos(theta); // x
      cities(i, 1) = R * sin(theta); // y
      outcities << setw(12) << i << setw(12) << cities(i, 0) << setw(12)
                << cities(i, 1) << endl;
    }
  }

  if(type == 1){

  for(int i=0;i<N;i++){

    cities(i,0)=_rnd.Rannyu(); // x
    cities(i,1)=_rnd.Rannyu(); // y
  outcities << setw(12) << i << setw(12) << cities(i, 0) << setw(12)
                << cities(i, 1) << endl;
  }

  }


  outcities.close();
  return cities;
}
mat Distances(const mat &cities) {
  // computes the distances between the generater cities
  mat D;
  int N = cities.n_rows;
  D.zeros(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // calcolo la distanza tra ogni elemento i,j
      D(i, j) = sqrt(pow(cities(i, 0) - cities(j, 0), 2) + pow(cities(i, 1) - cities(j, 1),2));
    }
  }
  return D;
}

///////////////////////////////////CHROMOSOME METHODS/////////////////////////////////////////////////////
Chromosome::Chromosome(int Ncities, Random &rnd) {

  _N_genes = Ncities;
  _genes.set_size(_N_genes);
  for (int i = 0; i < _N_genes; i++) {
    _genes(i) = i; // istanzio il primo come 1,2,3,4,5... e ottengo tutti gli
                   // altri da permutazioni (CC)
  }

  _rnd = rnd;
}

Chromosome::Chromosome() {
  _N_genes = 0;
  _genes.reset(); // empty constructor
}

bool Chromosome::Check_Boundary() {
  bool bound = true; // booleana per dire quando le cc sono soddisfatte

  // confronto tutti i geni nel cromosoma, se ne trovo due uguali, aumento il
  // contatore e do un altro valore random al gene uguale
  for (int i = 0; i < _N_genes; i++) {
    for (int j = 0; j < _N_genes; j++) {
      if (i != j && _genes(i) == _genes(j))
        bound = false;
    }
  }
  return bound;
}

void Chromosome::Swap(int Ntimes) {
  for (int i = 0; i < Ntimes; i++) {
    // select two random genes
    int index1 = (int)(1 + _rnd.Rannyu() * (_N_genes - 1));
    int index2 = (int)(1 + _rnd.Rannyu() * (_N_genes - 1));

    // swap them
    double a = _genes(index1);
    _genes(index1) = _genes(index2);
    _genes(index2) = a;
  }
  if (!this->Check_Boundary()) {
    cerr << "ERRORE: cromosoma invalido dopo Swap" << endl;
    exit(1);
  }
}

void Chromosome::Swap_genes(int i, int j) {
  double temp;
  temp = _genes(i);
  _genes(i) = _genes(j);
  _genes(j) = temp;
}

void Chromosome::Set_Genes(const vec &genes) {
  if (genes.n_elem != _N_genes) {
    cerr << "Errore: Set_Genes dimensione non corrispondente (" << genes.n_elem
         << " vs " << _N_genes << ")" << endl;
    exit(1);
  }

  for (int i = 0; i < _N_genes; i++) {
    if (genes(i) >= _N_genes) {
      cerr << "ERRORE: genes(" << i << ") = " << genes(i) << " >= " << _N_genes
           << endl;
      exit(1);
    }

    _genes(i) = genes(i);
  }
}

///////////////////////////////////GENETIC ALGORITHM METHODS///////////////////////////////////////

GA::GA(int N_chromosomes, int N_cities) {
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
        _rnd.SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;
  _rnd.SaveSeed();

  _N_chromosomes = N_chromosomes;

  _population.resize(_N_chromosomes);

  for (int i = 0; i < N_chromosomes; i++) {
    Chromosome chromo(N_cities, _rnd); // template chromosome
    chromo.Swap((int)N_cities);        // N swaps each time
    if (!chromo.Check_Boundary()) {    // Sempre valido
      cerr << "Errore: cromosoma iniziale invalido" << endl;
      exit(1);
    }
    _population[i] = chromo;
  }
  // set mutation and crossover probabilities
  _P_pair = 0.12;
  _P_shift = 0.2;
  _P_perm = 0.12;
  _P_invert = 0.15;
  _P_crossover = 0.5;
}

void GA::Swap_in_population(int i, int j) {
std::swap(_population[i], _population[j]);
}

int GA::GetN() { return _N_chromosomes; }

double GA::Loss(const mat &D, int i) { // prende in input il matricione che contiene ogni
                         // distanza (tanto sono fix), che calcolo nel main
  double len{};

  // salvo il percorso e sommo le distanze
  const Chromosome &c = _population[i];
  const vec &path = c.Get_Genes(); //[1,4,3...]
  int stops = c.Get_Ngenes() - 1;  // numero di fermate
  for (int j = 0; j < stops; j++) {
    int a = path(j);
    int b = path(j + 1);
    if (a < 0 || a >= D.n_rows || b < 0 || b >= D.n_cols) {
      cerr << "ERRORE: path(" << j << ")=" << a << ", path(" << j + 1
           << ")=" << b << " fuori dai limiti di D" << endl;
      exit(1);
    }
    // sommo tutte le distanze fino alla penultima
    len += D(a, b);
  }
  // aggiungo la distanza per tornare alla prima città
  len += D(path(stops), path(0));

  return len;
}

// ordina i cromosomi nella popolazione da quello con la loss minore a quello
// con la loss maggiore
void GA::Sort_by_fitness(const mat &D) {

  int N = GetN();
  std::vector<double> losses(N); //precalculate all losses once
  for (int i = 0; i < N; ++i)
    losses[i] = Loss(D, i);
//selection sort
  for (int i=0;i<N-1; i++) {
    int best=i;
    for (int j=i+1; j<N; j++) {
      if (losses[j] < losses[best]) {
        best=j;
      }
    }
    if (best != i) {
      std::swap(losses[i], losses[best]);      
      Swap_in_population(i, best);            
    }
  }

}

int GA::Selection(const mat &D) { // unfair roulette, picks the best paths with
                                  // greater probability. since the best ones
                                  // are in the lowest positions, i choose p=2

  int j = (int)(pow(_rnd.Rannyu(), 2) * this->GetN());
  return j; // returns the index corresponding to the selected chromosome
}

void GA::Pair_swap(Chromosome &c) { c.Swap(1); }

void GA::Shift_contiguous(Chromosome &c) { //shifts a random chunk of genes in the chromosome by a random amount of positions
  int N = c.Get_Ngenes();
  int start=(int)(_rnd.Rannyu(0.0, double(N-2)));      
  int length=(int)(_rnd.Rannyu(1,(N-start)/2.0)); 

  vec chunk(length); //store chunk of genes to be moved
  for(int i=0;i<length;i++){
    chunk(i)=c.Get_Gene(i+start);
  }
  //replace with the contiguous chunk
  for(int i=0;i<length;i++){
    c.Set_Gene(i+start,c.Get_Gene(i+start+length));
  }
  //overwrite the contiguous chunk
  for(int i=0;i<length;i++){
    c.Set_Gene(i+start+length, chunk(i));
  }

  if (!c.Check_Boundary()) {
    cerr << "ERRORE: cromosoma invalido dopo shift contiguo!" << endl;
    exit(1);
  }
}

void GA::Permutation_contiguous(Chromosome &c) {

  int N = c.Get_Ngenes();
  int begin = 1 + (int)_rnd.Rannyu() * (N - 1); // can't touch the first city
  int Nperm = (int)(_rnd.Rannyu() * (N - begin));
  for (int j = begin; j < Nperm + begin; j++) {
    c.Swap_genes(j, N - j - 1);
  }
  if (!c.Check_Boundary()) {
    cerr << "ERRORE: cromosoma invalido dopo permutation!" << endl;
    exit(1);
  }
}

void GA::Invert_path(Chromosome &c) {
  // genero le posizioni da swappare
  int N = c.Get_Ngenes();
  int m= (int)(_rnd.Rannyu(1.0,(double)(N-2))); //choose the number of contiguous cities to be inverted in the path
  int start = (int)_rnd.Rannyu(1.0,(double)(N-m)); // choose a random starting city

 
  for (int i=0; i<(int)(m/2);i++) {
    c.Swap_genes(i+start, i+start+m-1-i);
  }
  if (!c.Check_Boundary()) {
    cerr << "ERRORE: cromosoma invalido dopo invert!" << endl;
    exit(1);
  }
}

void GA::Random_mutation(
    Chromosome &c) { // selects a random mutation and performs it

  if (_rnd.Rannyu() <= _P_pair)
    Pair_swap(c);

 else if (_rnd.Rannyu() <= _P_shift)
    Shift_contiguous(c);

 else if (_rnd.Rannyu() <= _P_perm)
    Permutation_contiguous(c);

 else if (_rnd.Rannyu() <= _P_invert)
    Invert_path(c);
}

void GA::Crossover(Chromosome &f1, Chromosome &f2) { 
  int N = f1.Get_Ngenes();
  int cut = floor(_rnd.Rannyu(1,N)); // choose a random cutting position

  vec father(N), mother(N);
//fills mom and dad
  for (int k = 0; k < N; k++) {
    father(k) = f1.Get_Gene(k);
    mother(k) = f2.Get_Gene(k);
  }

  vec son1(N), son2(N);
  //fills sons' chromosomes until cutting position is reached
  for (int k = 0; k < cut; k++) {
    son1(k) = father(k);
    son2(k) = mother(k);
  }
//fills son 1 with mother's genes without repetitions
  int pos1 = cut;
  for (int k = 0; k < N; k++) {
    bool already_in = false;
    for (int m = 0; m < pos1; m++) {
      if (mother(k) == son1(m)) {
        already_in = true;
        break;
      }
    }
    if (!already_in) {
      son1(pos1) = mother(k);
      pos1++;
    }
  }

  int pos2 = cut;
  for (int k = 0; k < N; k++) {
    bool already_in = false;
    for (int m = 0; m < pos2; m++) {
      if (father(k) == son2(m)) {
        already_in = true;
        break;
      }
    }
    if (!already_in) {
      son2(pos2) = father(k);
      pos2++;
    }
  }

  f1.Set_Genes(son1);
  f2.Set_Genes(son2);
  if (!f1.Check_Boundary() || !f2.Check_Boundary()) {
    cerr << "error! the algorithm generated a chromosome that doesn't respect "
            "boundary conditions. Check the code"
         << endl;
    for (int i = 0; i < 34; i++) {
      cerr << "\t" << f1.Get_Gene(i);
    }
    exit(1);
  }
}

void GA::Evolve(const mat &D) {

  // sort the population by fitness
  Sort_by_fitness(D);

  // create a new population
  vector<Chromosome> new_population;
  //elitism
  int elitism=1;
  for(int i=0;i<elitism;i++)
  new_population.push_back(_population[i]);
  // loop over current population's chromosomes, while leaving one spot free for elitism
  while (new_population.size() < _N_chromosomes) {
    // choose two parents
    int mom = Selection(D);
    int dad = Selection(D);
    // generetes two sons
    Chromosome f1 = _population[mom];
    Chromosome f2 = _population[dad];
    // crossover
    if (_rnd.Rannyu() <= _P_crossover)
      Crossover(f1, f2);

    // random mutations
    Random_mutation(f1);
    Random_mutation(f2);

    new_population.push_back(f1);
    new_population.push_back(f2);
  }
// Se manca un solo cromosoma per arrivare a N, ne aggiungo uno solo
  if ((int)new_population.size() < _N_chromosomes) {
    int parent = Selection(D);
    Chromosome f = _population[parent];
    Random_mutation(f);
    new_population.push_back(f);
  }
  // replace the old population with the new one
  _population = new_population;
}

void GA::Pick_best(const mat &D) {
  Sort_by_fitness(D);
  vec best = _population[0].Get_Genes();
  ofstream out;
  out.open("../OUTPUT/best_path.dat");
  for (int i = 0; i < best.size(); i++) {
    out << best(i) << endl;
  }
  out.close();

  cout<<"best length="<<Loss(D,0);
}


double GA::BestHalfLoss(const mat &D){

double acc=0;
 // Sort_by_fitness(D);
  for(int i=0;i<(int)GetN()/2;i++){
    acc+=Loss(D,i);
  }
  acc/=((int)GetN()/2);
  return acc;

}
