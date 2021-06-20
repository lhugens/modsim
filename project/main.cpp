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

void NVT_high_density(){

    simulation s("output/nvt_high_density/coords", 4.0);
    s.dl = 0.5;
    s.dn = 0.5;
    s.fcc_config(5, 0.65);
    s.print_parameters();
    s.write_config(0);
    assert(!(s.exists_general_overlap()));

    int total_steps = 1e6;
    int fast_steps = 1e3;
    int measure_steps = 1e3;

    while(s.step < total_steps){
        s.NVT_run(fast_steps);
        s.update_S();
        s.write_config(s.step);
        cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%]" 
            << " acc. rate: " << setw(10) << s.accept_rate 
            << " rho: " << setw(10) << s.rho
            << " S: " << setw(10) << s.S
            << flush;
    }
    cout << endl;
}

void S_L_plot(){
    double L_min = 3.0;
    double L_max = 5.0;
    int L_no  = 5;
    double delta = (L_max - L_min) / (double)L_no;

    int total_steps = 1e6;
    int fast_steps = 1e3;
    int measure_steps = 1e3;

    //for(int i=0; i<L_no; i++){
    int i = 1;
    string no_str = to_string(i);
    string filename = "output/SL/SL_" + string(3 - no_str.length(), '0') + no_str + ".dat";
    string folder = "output/coords_SL/coords" + string(3 - no_str.length(), '0') + no_str;
    ofstream file;
    file.open(filename);

    double L = L_min+delta*i;
    L = 50;

    file << L << endl;

    simulation s(folder, L);
    s.dl = 0.1;
    s.dn = 0.5;
    s.fcc_config(5, 1.0);
    s.print_parameters();
    s.write_config(0);
    assert(!(s.exists_general_overlap()));

    while(s.step < total_steps){
        s.NVT_run(fast_steps);
        s.update_S();
        s.write_config(s.step);
        cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%]" 
            << " acc. rate: " << setw(10) << s.accept_rate 
            << " S: " << setw(10) << s.S
            << flush;
    }
    cout << endl;

    for(int i=0; i<measure_steps; i++){
        s.step++;
        s.NVT_step();
        s.update_S();
        file << left << setw(20) << i << setw(20) << s.S << endl;
    }
    file.close();

    //}
    
}


void EOS(){
    // preparing betaPv0 values
    double betaPv0_min = 2.0;
    double betaPv0_max = 12.0;
    int betaPv0_no  = 10;
    double delta = (betaPv0_max - betaPv0_min) / (double)betaPv0_no;

    // steps control
    int total_steps = 1e6;
    int fast_steps = 1e3;
    int measure_steps = 1e4;

    simulation s("", 3);
    s.file_config("liquid_0.374.dat");

    // MC proposals
    s.dl = 0.01;
    s.dn = 0.1;
    s.dV = 10;
    s.pvol = 0.1;

    for(int i=0; i<betaPv0_no; i++){

        // prepare file to output density and S values
        string no_str = to_string(i);
        string filename = "EOS/EOS_" + string(3 - no_str.length(), '0') + no_str + ".dat";

        ofstream file;
        file.open(filename);
        double betaPv0 = betaPv0_min+delta*i;
        file << betaPv0 << endl;

        // start a new simulation, using the final config of
        // the previous simulation as initial config
        s.step = 0;
        s.accept = 0;
        s.folder = "coords_EOS/coords_" + string(3 - no_str.length(), '0') + no_str;

        s.betaP = betaPv0 / s.v0;
        s.print_parameters();
        s.write_config(0);

        while(s.step < total_steps){
            s.NPT_run(fast_steps);
            s.write_config(s.step);

            cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%]" 
                 << " acc. rate: " << setw(10) << s.accept_rate
                 << " rho: "       << setw(10) << s.rho
                 << flush;
        }
        cout << endl;

        for(int i=0; i<measure_steps; i++){
            s.step++;
            s.NPT_step();
            s.update_S();
            file << left << setw(20) << i << setw(20) << s.rho << setw(20) << s.S << endl;
        }
        file.close();
    }
}

void NPT_1_value(double betaPv0, string f){

    string folder = "output/" + f + "/coords";
    simulation s(folder, 3);
    s.file_config("configs/liquid_0.374.dat");
    s.betaP = betaPv0 / s.v0;

    string filename = "output/" + f + "/eos";
    ofstream file;
    file.open(filename);
    file << betaPv0 << endl;

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

void homogenize(){
    simulation s("output/homogenize/coords", 1.0);
    s.dl = 0.5;
    s.dn = 0.5;
    s.file_config("configs/liquid_0.374.dat");
    s.print_parameters();
    assert(!(s.exists_general_overlap()));

    int total_steps = 1e6;
    int fast_steps = 1e3;
    int measure_steps = 1e3;

    while(s.step < total_steps){
        s.NVT_run(fast_steps);
        s.update_S();
        s.write_config(s.step);
        cout << "\r [" << setw(3) << round((double)s.step * 100 /total_steps)  << "%]" 
            << " acc. rate: " << setw(10) << s.accept_rate 
            << " rho: " << setw(10) << s.rho
            << " S: " << setw(10) << s.S
            << flush;
    }
    cout << endl;
}


int main(){
    dsfmt_seed(time(NULL));
    //EOS();
    //NPT_1_value(0.25, "NPT4");
    S_L_plot();
    //homogenize();
    //NVT_high_density();
    // 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
}
