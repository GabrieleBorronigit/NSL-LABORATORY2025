#pragma once
#include <cmath>
#include <fstream>
#include <iostream>
#include "random.h"
#include <string>
using namespace std;



//create template function class, useful for scalability
class Template_function {
 public:
  // since this is an abstract class, it doesn't need any constructor because I never call an object of this class explicitly
  virtual double Eval(
      double x) const = 0;  // Pure-virtual method for f(x)

  virtual ~Template_function() {};
};



// requested integrand function
class Integrand : public Template_function {
 public:
  double Eval(double x) const override {
    return M_PI / 2.0 * cos(x * M_PI / 2.0);
  }
};





// Class for probability distribution functions, it has an extra method for sampling with the inverse of the cumulative of p(x)
class PDF : public Template_function {
 public:
  virtual double Sample(Random &gen) const = 0;
};

// samples as p(x)=Integrand
class IS_P_integrand : public PDF {
 public:
  double Eval(double x) const override {
    return M_PI / 2.0 * cos(x * M_PI / 2.0);
  }
  double Sample(Random &gen) const override {
    return 2 / M_PI * asin(gen.Rannyu());
  }
};

// Linear pdf sampling: p(x)=2(1-x),  obtained by observing that integrand is a decreasing function in [0,1] 
//that shows a maximum in x=0 and it is null in x=1
// so P(x)=A(1-x) and normalize to find A=2

class IS_P_linear : public PDF {
 public:
  double Eval(double x) const override { return 2 * (1 - x); }
  double Sample(Random &gen) const override {
  
    return 1 - sqrt(1 - gen.Rannyu());
  }
};




//class for MC integrals
class Integral_MC {



 private:                                       
  double m_a, m_b;                              //integration interval
  int m_N;                                      //number of samples
  Random m_gen;                                 //RNG



 public:

  Integral_MC(double min, double max, int punti, Random rnd) //constructor
      : m_a{min}, m_b{max}, m_N{punti}, m_gen{rnd} {} ;




  // metodi per accedere agli attributi
  void Set_interval(double min, double max) {
    m_a = min;
    m_b = max;
  }
  void Set_points(int npoints) { m_N = npoints; }

  double GetA() { return m_a; }
  double GetB() { return m_b; }
  double GetN() { return m_N; }





  // metodo per integrale MC con distribuzione uniforme (metodo della media)
  double Uniform_MC(const Template_function &f);



  // metodo per integrale MC con importance sampling
  // per generalizzare un minimo, il metodo accetta l'integranda e la funzione
  // per campionare i punti secondo la d.d.p. voluta
  double Importance_Sampling(const Template_function &f, PDF &p);
};




//funzioni per analisi dati

double error(double av, double av2,
             int n);  // ritorna la deviazione standard della media

// funzione che calcola media e deviazione standard usando data blocking

void ComputeMeanErrors_DB(const vector<double> &data, int N_blocks,
                          string filename);