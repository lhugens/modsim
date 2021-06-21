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

void NPT(double betaPv0){

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

void NPT_continue(double betaPv0){
    string b = to_string(betaPv0);
    string folder = "output/NPT/NPT_" + string(10 - b.length(), '0') + b;
    fs::current_path(folder);
    string filename, last_file, counter_str;

    cout << fs::current_path() << endl;

    for (const auto& dirEntry : fs::recursive_directory_iterator(fs::current_path())){
        filename = dirEntry.path().string();

        if(filename.find("coords") != string::npos){
            last_file = filename;
        }
    }
    cout << last_file << endl;

    unsigned first = last_file.find("p_");
    unsigned last = last_file.find(".dat");
    counter_str = last_file.substr(first+2, last-first-2);

    int counter = stoi(counter_str);

    cout << counter << endl;


    ofstream file;
    file.open("eos.dat", std::ios_base::app);

    
    // initialize simulation
    simulation s("coords/coords", 3);
    s.file_config(last_file);
    s.betaP = betaPv0 / s.v0;

    s.dl = 0.01;
    s.dn = 0.1;
    s.dV = 10;
    s.pvol = 0.1;
    s.print_parameters();

    int total_steps = 1e6;
    int fast_steps = 1e4;

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
    file.close();
}


/*
void convert(){
    cout << filesystem::current_path() << endl;
    //using recursive_directory_iterator = filesystem::recursive_directory_iterator;
    string filename;

    for (const auto& dirEntry : filesystem::recursive_directory_iterator(filesystem::current_path())){
        filename = dirEntry.path().string();
        if(filename.find(".dat") != string::npos){
            cout << "converting " << filename  << endl;
            auto [part, box_l, N, L] = file_config_initial_old(filename);
            write_config_to_file(part, box_l, N, L, filename);
        }
    }
}
*/

int main(int argc, char **argv){
    dsfmt_seed(time(NULL));

    // getting betaPv0 value from command line argument
    double betaPv0 = stod(argv[1]);

    NPT_continue(betaPv0);


}







