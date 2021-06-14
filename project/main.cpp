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

void NVT(double rho, double L){
    // create simulation and fix length L
    simulation s(L);

    // create initial fcc configuration with (N_side, rho)
    s.fcc_config(5, rho);

    // write to config to file
    s.write_config(1);

    // make sure there's no overlap at the start
    assert(!(s.exists_overlap()));

    // main monte carlo loop
    int total_steps = 1000; 
    int equilibration_steps = (int)((double)total_steps/5);

    for(s.step=1; s.step<equilibration_steps; s.step++){
        s.propose_NVT();
        s.metropolis_acceptance_NVT();
        s.write_config(s.step);
        cout << "acceptance rate: " << s.accept_rate << endl;
    }
}

void NPT(double betaP, double L){
    // create simulation and fix length L
    simulation s(L);
    s.betaP = betaP;

    // create initial fcc configuration with (N_side, rho)
    s.fcc_config(5, 0.1);

    // write to config to file
    s.write_config(1);

    // make sure there's no overlap at the start
    assert(!(s.exists_overlap()));

    // main monte carlo loop
    int total_steps = 1000; 
    int equilibration_steps = (int)((double)total_steps/5);

    for(s.step=1; s.step<equilibration_steps; s.step++){
        s.propose_NPT();
        s.metropolis_acceptance_NPT();
        s.write_config(s.step);
        cout << "acceptance rate: " << s.accept_rate << endl;
    }
}

int main(){
    dsfmt_seed(time(NULL));
    NVT(0.5, 3.5);

}
