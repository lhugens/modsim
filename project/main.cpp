#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>

#define ndim 3

using namespace std;

#include "mt19937.h"
#include "tools.hpp"
#include "simulation.hpp"

int main(){
    dsfmt_seed(time(NULL));

    // create simulation and fix length L
    simulation s(3.5);

    // create initial fcc configuration with (N_side, rho)
    s.fcc_config(5, 1);

    // write to config to file
    s.write_config(1);

    // make sure there's no overlap at the start
    assert(!(s.exists_overlap()));

    /*
    // main monte carlo loop
    int total_steps = 1000; 
    int equilibration_steps = (int)((double)total_steps/5);

    for(int step=1; step<equilibration_steps; step++){
        s.propose_dl_dn();
        s.metropolis_acceptance();
        s.write_config(step);
        cout << "acceptance rate: " << s.accept_rate << endl;
    }

    for(int step=equilibration_steps; step<total_steps; step++){
        s.propose_dl_dn();
        s.metropolis_acceptance();
        s.write_config(step);
        cout << "acceptance rate: " << s.accept_rate << endl;
    }
    */

}
