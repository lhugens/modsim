#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <string>
#include <filesystem>
#include <string.h>

#define ndim 3

using namespace std;

namespace fs = filesystem;

#include "mt19937.h"
#include "tools.hpp"
#include "simulation.hpp"

int main(int argc, char **argv){
    dsfmt_seed(time(NULL));
    
    // getting betaPv0 value from command line argument
    double betaPv0 = stod(argv[1]);

    // preparing folder and file to store output
    string b = to_string(betaPv0);
    string folder = "output/NPT/NPT_" + string(10 - b.length(), '0') + b;
    string coords_folder = folder + "/coords";
    string coords_file   = coords_folder + "/coords";
    string filename = folder + "/eos.dat";
    fs::create_directory(folder);
    fs::create_directory(coords_folder);
    ofstream file;
    file.open(filename);
    file << betaPv0 << endl;

    // initialize simulation
    simulation s(coords_file, 3);
    s.file_config("configs/liquid_0.374.dat");
    s.betaP = betaPv0 / s.v0;

    s.dl = 0.01;
    s.dn = 0.1;
    s.dV = 10;
    s.pvol = 0.1;
    s.print_parameters();
    s.write_config(0);

    int total_steps = 1e6;
    int fast_steps = 1e4;
    int counter = 0;

    while(true){
        counter++;

        s.NPT_run(fast_steps);
        s.write_config(counter);
        s.update_S();

        file << left << setw(20) << counter << setw(20) << s.rho << setw(20) << s.S << endl;

        cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%]" 
             << " acc. rate: " << setw(10) << s.accept_rate
             << " rho: "       << setw(10) << s.rho
             << flush;
    }
    cout << endl;
}














