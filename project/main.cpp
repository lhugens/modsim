#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <string>

#define ndim 3

using namespace std;

#include "mt19937.h"
#include "tools.hpp"
#include "simulation.hpp"

void NVT(string folder, double rho, double L, int total_steps, int write_steps){
    // create simulation and fix length L
    simulation s(folder, L);

    // create initial fcc configuration with (N_side, rho)
    s.fcc_config(5, rho);

    // write to config to file
    s.write_config(0);

    // make sure there's no overlap at the start
    assert(!(s.exists_general_overlap()));

    while(s.step < total_steps){
        for(int i=0; i<write_steps; i++){
            s.step++;
            s.propose_NVT();
            s.metropolis_acceptance_NVT();
        }
        s.write_config(s.step);
        //cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps) << "%] " << "general: " << s.exists_general_overlap() << flush;
        cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%]" << " acc. rate: " << s.accept_rate 
                                                                                       << flush;
    }
    cout << endl;

}

void NPT(double betaP, double L){

    /*
    int total_steps = 1e4; 
    int measuring_steps = 100;
    int equilibration_steps = total_steps - measuring_steps;

    for(s.step=1; s.step<equilibration_steps; s.step++){
        s.propose_NVT();
        s.metropolis_acceptance_NVT();
        cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%] " << "acc. rate: " << s.accept_rate << flush;
    }
    cout << endl;
    cout << " WRITING TO FILE !!! " << endl;

    for(s.step=equilibration_steps; s.step<total_steps; s.step++){
        s.propose_NVT();
        s.metropolis_acceptance_NVT();
        s.write_config(s.step);
    }
    */
}

int main(){
    dsfmt_seed(time(NULL));
    NVT("coords_NEMATIC", 0.5, 4.8, 1e7, 1e4);
}
