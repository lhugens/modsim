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
    assert(!(s.exists_initial_overlap()));

    // main monte carlo loop
    int total_steps = 1e4; 
    int measuring_steps = (int)((double)total_steps/5);
    measuring_steps = 100;
    int equilibration_steps = total_steps - measuring_steps;

    for(s.step=1; s.step<equilibration_steps; s.step++){
        s.propose_NVT();
        s.metropolis_acceptance_NVT();
        cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%] " << "acc. rate: " << s.accept_rate << flush;
    }
    cout << endl;
    cout << " WRITING TO FILE !!! " << endl;
    cout << " WRITING TO FILE !!! " << endl;
    cout << " WRITING TO FILE !!! " << endl;
    cout << " WRITING TO FILE !!! " << endl;
    cout << " WRITING TO FILE !!! " << endl;

    for(s.step=equilibration_steps; s.step<total_steps; s.step++){
        s.propose_NVT();
        s.metropolis_acceptance_NVT();
        s.write_config(s.step);
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
    assert(!(s.exists_initial_overlap()));

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
