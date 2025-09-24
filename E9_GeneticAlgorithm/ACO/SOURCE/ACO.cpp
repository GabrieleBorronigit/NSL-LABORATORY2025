
#include "ACO.h"
#include "random.h"
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;
using namespace arma;

////////////////////////INITIALIZATION OF TSP METHODS///////////////////////////
mat Load_cities(string filename, int Ncities){

  ifstream in(filename);
  if(!in.is_open()){
    cerr<<"cannot open file "<<filename<<endl;
    exit(1);
  }
  mat cities(Ncities, 2);
  cities.zeros(Ncities,2);
    for(int i=0;i<Ncities;i++){
      double x{}, y{};
      in>>x;
      in>>y;
      cities(i,0)=x;
      cities(i,1)=y;
    }
    in.close();
    return cities;
}
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

  if (type == 1) {

    for (int i = 0; i < N; i++) {

      cities(i, 0) = _rnd.Rannyu(); // x
      cities(i, 1) = _rnd.Rannyu(); // y
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
      D(i, j) = sqrt(pow(cities(i, 0) - cities(j, 0), 2) +
                     pow(cities(i, 1) - cities(j, 1), 2));
    }
  }
  return D;
}

////////////////////ANT METHODS/////////////////////////////

Ant::Ant(int Ncities, Random &rnd) {
  Ncities_ = Ncities;
  visited_.push_back(0);
  for (int i = 1; i < Ncities; i++) {
    unvisited_.push_back(i);
  }
  rnd_ = rnd;
}



void Ant::Walk(const mat &Feromons, const mat &D, double alpha, double beta) {

  Pij_.clear();
  // calcolo prob di andare in ogni cita non ancora vista
  int current = visited_.back(); // ultima citta visitata
  // calcolo tutte le probabilità
  for (int j = 0; j < unvisited_.size(); j++) {
    int next = unvisited_[j];
    Pij_.push_back(pow(1 / D(current, next), alpha) *
                   pow(Feromons(current, next), beta));
  }
  // normalizzo
  double norm = accumulate(Pij_.begin(), Pij_.end(), 0.0);
  for (int j = 0; j < Pij_.size(); j++) {
    Pij_[j] /= norm;
  }
  // choose next city index
  int step = Next_step();

  // perform the step
  visited_.push_back(unvisited_[step]);
  unvisited_.erase(unvisited_.begin() + step);
}

int Ant::Next_step() { // sceglie il prox step con probabilità pesata dai
                       // feromoni
  double r = rnd_.Rannyu();
  double cumulative = 0.0;

  for (int i = 0; i < Pij_.size(); i++) {
    cumulative += Pij_[i];

    if (r <= cumulative) // cumulativa(j-1)<j<cumulativa(j)
      return i;
  }
  // fallback di sicurezza, nel caso di arrotondamenti
  return Pij_.size() - 1;
}

void Ant::Reset() {
  visited_.clear();
  unvisited_.clear();
  Pij_.clear();

  visited_.push_back(0); // restart from fisrst city
  for (int i = 1; i < Ncities_; i++)
    unvisited_.push_back(i);
}

vector<int> Ant::Get_visited() { return visited_; }

int Ant::Get_visited(int i) { return visited_[i]; }

/////////////////////////COLONY METHODS/////////////////////////////

Colony::Colony(int Ncities, int Nants, double rho, double Q) {

  // initialize random numbers generator
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
        rnd_.SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;
  rnd_.SaveSeed();
  // colony parameters
  Ncities_ = Ncities;
  Nants_ = Nants;
  // initialize colony
  for (int i = 0; i < Nants_; i++) {
       int s[4];
    for (int j = 0; j < 4; j++) s[j] = int(1e4 * rnd_.Rannyu());

    Random r;
    r.SetRandom(s, p1, p2);
    Ant ant(Ncities_, r);
    Colony_.push_back(ant);
  }
  // initialize feromones (for the first step it is equal for every path)
  Feromons_.set_size(Ncities_, Ncities_);
  Feromons_.fill(1.0);
  Feromons_.diag().zeros();  // diagonal elements set to 0 (can t stay in the same city)


  // evaporation rate
  rho_ = rho;
  // feromones deposition rate
  Q_ = Q;
}

double Colony::Loss(const mat &D, int i) {

  double len{};
  // salvo il percorso e sommo le distanze
  const vector<int> &path = Colony_[i].Get_visited(); //[1,4,3...]
  int stops = Ncities_ - 1;                           // numero di fermate
  for (int j = 0; j < stops; j++) {
    int a = path[j];
    int b = path[j + 1];
    if (a < 0 || a >= D.n_rows || b < 0 || b >= D.n_cols) {
      cerr << "ERRORE: path(" << j << ")=" << a << ", path(" << j + 1
           << ")=" << b << " fuori dai limiti di D" << endl;
      exit(1);
    }
    // sommo tutte le distanze fino alla penultima
    len += D(a, b);
  }
  // aggiungo la distanza per tornare alla prima città
  len += D(path[stops], path[0]);

  return len;
}

void Colony::Update_feromons(const mat &D) {

  // evaporate previous feromons
  for (int i = 0; i < Ncities_; i++) {

    for (int j = 0; j < Ncities_; j++) {
      Feromons_(i, j) *= (1 - rho_);
    }
  }

  // deposit the new ones tauij=tauij+sum{Dtauij}
  Sort_by_fitness(D); //elitism: the top 40% of ants guides the evolution
  for (int a = 0; a < (int)(Nants_*0.4); a++) {
    const vector<int> &path = Colony_[a].Get_visited();
    double loss = Loss(D, a);
    for (int i = 0; i < Ncities_ - 1; i++) {
      // Update the feromons only on the corresponding arc!
      int arc_start = path[i];
      int arc_end = path[i + 1];
      Feromons_(arc_start, arc_end) += Q_ / loss;
    }
    // close the trip by returning to the first city
    Feromons_[path.back(), path.front()] += Q_ / loss;
  }
  for (int i = 0; i < Ncities_; i++) {
    Feromons_(i, i) = 0.0;
  }
}

void Colony::Search(const mat &D, double alpha, double beta) {
  // let the ants walk
  for (int a = 0; a < Nants_; a++) {
    Colony_[a].Reset();
    for (int i = 0; i < Ncities_; i++) {
      Colony_[a].Walk(Feromons_, D, alpha, beta);
    }
  }
  // deposit feromons to drive the next walkers
  Update_feromons(D);
}

void Colony::Sort_by_fitness(const mat &D) {

  int N = Nants_;
  std::vector<double> losses(N); // precalculate all losses once
  for (int i = 0; i < N; ++i)
    losses[i] = Loss(D, i);
  // selection sort
  for (int i = 0; i < N - 1; i++) {
    int best = i;
    for (int j = i + 1; j < N; j++) {
      if (losses[j] < losses[best]) {
        best = j;
      }
    }
    if (best != i) {
      std::swap(losses[i], losses[best]);
      Swap_in_colony(i, best);
    }
  }
}

void Colony::Swap_in_colony(int i, int j) { std::swap(Colony_[i], Colony_[j]); }

void Colony::Pick_best(const mat &D) {
  ofstream out("../OUTPUT/best_path.dat");
  Sort_by_fitness(D);
  for (int i = 0; i < Ncities_; i++) {
    out << Colony_[0].Get_visited(i) << endl;
  }
  out.close();
  cout << "best length=" << Loss(D, 0)<<endl;
}

double Colony::Best_loss(const mat &D) {
  Sort_by_fitness(D);
  return Loss(D, 0);
}

double Colony::BestHalfLoss(const mat &D) {
  // Sort_by_fitness(D);
  double acc{};
  for (int i = 0; i < (int)Nants_ / 2; i++) {
    acc += Loss(D, i);
  }
  acc /= ((int)Nants_ / 2);
  return acc;
}
