#include <iostream>
#include "system.h"

using namespace std;

int main (){
  int nconf = 1;
  System Sys;
  Sys.initialize();
  Sys.initialize_properties();
  Sys.block_reset(0);

  for(int i=0; i < Sys.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < Sys.get_nsteps(); j++){ //loop over steps in a block
      Sys.step();
      Sys.measure();
      if(j%50 == 0){
        Sys.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++;
      }
    }
    Sys.averages(i+1);
    Sys.block_reset(i+1);
  }

  Sys.finalize();

  return 0;
}