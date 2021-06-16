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
        cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%]" << " acc. rate: " << s.accept_rate << flush;
    }
    cout << endl;

}

void NPT(int betaPv0_no, string folder, double betaPv0, double rho_initial, double L, int total_steps, int write_steps, int measure_steps){

    string no_str = to_string(betaPv0_no);
    string filename = "EOS/EOS_" + string(3 - no_str.length(), '0') + no_str + ".dat";

    ofstream file;
    file.open(filename);
    file << betaPv0 << endl;

    // create simulation and fix length L
    simulation s(folder, L);
    s.fcc_config(5, rho_initial);
    s.betaP = betaPv0 / s.v0;
    s.write_config(0);
    cout << "betaP " << s.betaP << endl;
    cout << "rho   " << s.rho << endl;

    // make sure there's no overlap at the start
    assert(!(s.exists_general_overlap()));

    while(s.step < total_steps){
        for(int i=0; i<write_steps; i++){
            s.step++;
            s.NPT_step();
        }
        s.write_config(s.step);

        // dont forget to take this out
        file << left << setw(20) << s.step << setw(20) << s.rho << endl;

        cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%]" 
             << " betaPv0: "   << setw(3) << betaPv0 
             << " acc. rate: " << setw(10) << s.accept_rate
             << " rho: "       << setw(10) << s.rho
             << flush;
    }
    cout << endl;

    for(int i=0; i<measure_steps; i++){
        s.step++;
        s.NPT_step();
        file << left << setw(20) << i << setw(20) << s.rho << endl;
    }

    file.close();
}

void EOS(){
    double betaPv0_min = 2.0;
    double betaPv0_max = 12.0;
    int betaPv0_no  = 10;
    double delta = (betaPv0_max - betaPv0_min) / (double)betaPv0_no;

    for(int i=0; i<betaPv0_no+1; i++){
        NPT(i, "coords_EOS", betaPv0_min+delta*i, 0.3, 3, 1e5, 1e3, 1e2);
    }
}

int main(){
    dsfmt_seed(time(NULL));
    NPT(0, "coords_EOS", 2.0, 0.3, 3, 1e4, 10, 0);
}






