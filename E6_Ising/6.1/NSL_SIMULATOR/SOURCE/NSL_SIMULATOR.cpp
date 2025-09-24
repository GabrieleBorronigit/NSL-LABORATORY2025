/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "system.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);
/*
//usando restart
  // equilibrio
  int equilibrium = 20000;
  // ciclo sulle temperature
  double Tmin = 0., Tmax = 2.0;
  for (double T = Tmin; T <= Tmax; T += 0.025) {

    SYS.read_configuration();
    SYS.setTemp(T);
     for(int i=0;i<equilibrium;i++){
    SYS.step();
  }
    for(int i=0; i<SYS.get_nbl();i++){

      for(int j=0;j<SYS.get_nsteps();j++){
        SYS.step();
        SYS.measure();

      }
    SYS.averages(i+1);
    SYS.block_reset(i+1);
    }
    SYS.finalize();
    }
*/

  // equilibrio
  int equilibrium = SYS.get_nsteps()/2.0;
  // ciclo sulle temperature
  double Tmin = 0.5, Tmax = 2.0;
  for (double T = Tmin; T <= Tmax; T += 0.025) {
    SYS.setTemp(T);

    for (int i = 0; i < equilibrium; i++) {
      SYS.step();
    }
    SYS.Reset_acceptance();
    // evolvo e misuro
    for (int i = 0; i < SYS.get_nbl(); i++) {      // loop over blocks
      for (int j = 0; j < SYS.get_nsteps(); j++) { // loop over steps in a block
        SYS.step();
        SYS.measure();
        if (j % 50 == 0) {
          // SYS.write_XYZ(nconf); //Write actual configuration in XYZ format
          // //Commented to avoid "filesystem full"!
         nconf++;
        }
      }
      SYS.averages(i + 1);
      SYS.block_reset(i + 1);}
      SYS.reset_acc();
    
  }
  SYS.finalize();

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

// metti potenziale deltiforme, inizializzando le velocita in mood che all fine
// il momento risultante globale sia zero es particella 1 vx=3 particellla 2
// vx=-3 ecc, scegli le velocità in base al legame con la temperatura!!! nel
// file pofv.dat ,, stampa anche le medie di blocco!! alla fine devi fare un
// grafico con la distribuzione delle velocità, facendo vedere che partendo da
// una delt,a, questa va ad allargarsi blocco dopo blocco