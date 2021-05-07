#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "mt19937.h"

using namespace std;

#define ndim 3

struct dynamics{
    int npart;                  // number of particles
    int step;                   // current time step of the simulation
    double box_l;               // (cubic) box side length
    double dt;                  // time step
    vector<vector<double>> r0;  // array to store particle positions at t - dt
    vector<vector<double>> r1;  // array to store particle positions at t

    dynamics(int NPART, double DT, double BOX_L) :
        npart(NPART), dt(DT), box_l(BOX_L)
    {
        r0.resize(npart, vector<double>(ndim));
        r1.resize(npart, vector<double>(ndim));
        dsfmt_seed(time(NULL));
    }

    inline double rand(){
        return dsfmt_genrand();
    }

    void randomize_positions(){
        for(int i=0; i<npart; i++){
            for(int j=0; j<ndim; j++){
                r0[i][j] = box_l * rand();
            }
        }
    }

    void write_positions_to_file(){
        char filename[128];
        sprintf(filename, "coords/coords_step%07d.dat", step);

        ofstream file;
        file.open(filename);

        file << npart << endl;

        for(int i=0; i<ndim; i++){
            file << left << setw(10) << 0.0 
                         << setw(10) << box_l
                         << endl;
        }

        for(int i=0; i<npart; i++){
            file << left;
            for(int j=0; j<ndim; j++){
                file << setw(10) << r0[i][j];
            }
            file << setw(10) << 0.5 << endl;
        }
        file.close();
    }

};

int main(){
    dynamics md(10, 0.01, 10);
    md.randomize_positions();
    md.write_positions_to_file();
}

