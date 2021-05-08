#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "mt19937.h"

using namespace std;

#define ndim 3

struct dynamics{
    int npart;                  // number of particles
    int step = 0;               // current time step of the simulation
    double box_l;               // (cubic) box side length
    double dt;                  // time step
    double m = 1;               // particle mass
    double beyta = 0.5;         // beyta = 1 / (k_B * T)
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

    void write_positions_to_file(){
        char filename[128];
        sprintf(filename, "coords/coords_step%07d.dat", step);

        ofstream file;
        file.open(filename);

        file << npart << endl;

        for(int i=0; i<ndim; i++){
            file << left << setw(20) << 0.0 
                         << setw(20) << box_l
                         << endl;
        }

        for(int i=0; i<npart; i++){
            file << left;
            for(int j=0; j<ndim; j++){
                file << setw(20) << r0[i][j];
            }
            file << setw(20) << 0.5 << endl;
        }
        file.close();
    }

    void print(vector<vector<double>> &v){
        for(int i=0; i<v.size(); i++){
            cout << left;
            for(int j=0; j<v[0].size(); j++){
                cout << setw(20) << v[i][j];
            }
            cout << endl;
        }
    }

    void print(vector<double> &v){
        for(int i=0; i<v.size(); i++){
            cout << v[i] << endl;
        }
    }

    double square_norm(vector<double> &v){
        double sum = 0;
        for(int i=0; i<ndim; i++){
            sum += pow(v[i], 2);
        }
        return sum;
    }

    vector<double> force(vector<double> &ri, vector<double> &rj){
        vector<double> rij(ndim);
        vector<double> fij(ndim);
        for(int k=0; k<ndim; k++){
            rij[k] = ri[k]-rj[k];
        }
        double sq = square_norm(rij);
        double factor = 48 * (1/pow(sq, 3) - 0.5) / pow(sq, 4);
        for(int k=0; k<ndim; k++){
            fij[k] = rij[k] * factor;
        }
        return fij;
    }


    void initialize_positions_velocities(){

        // random r0 positions, each component in [-box_l, box_l]
        for(int i=0; i<npart; i++){
            for(int j=0; j<ndim; j++){
                r0[i][j] = box_l * rand();
            }
        }

        // randomize velocities, each component in [-1,1]
        // and compute total momentum
        vector<vector<double>> v0(npart, vector<double>(ndim));
        vector<double> P(ndim);

        for(int j=0; j<ndim; j++){
            for(int i=0; i<npart; i++){
                v0[i][j] = rand() * 2 - 1;
                P[j] += v0[i][j];
            }
        }

        // shift velocities s.t. total momentum is 0
        for(int j=0; j<ndim; j++){
            for(int i=0; i<npart; i++){
                v0[i][j] -= P[j] / (double)npart;
            }
        }

        // scale velocities s.t. E_kinetic = 3 N k_B T / 2
        double sum_square_speeds = 0;
        for(int i=0; i<npart; i++){
            sum_square_speeds += square_norm(v0[i]);
        }
        double scaling_factor = sqrt(3 * npart / (beyta * m * sum_square_speeds));

        for(int i=0; i<npart; i++){
            for(int j=0; j<ndim; j++){
                v0[i][j] *= scaling_factor;
            }
        }

        // determine r1 based on these velocities 
        // and on calculated forces
        double force_coefficient = pow(dt, 2) / (2*m);
        vector<double> fij;

        for(int i=0; i<npart-1; i++){
            for(int j=i+1; j<npart; j++){
                fij = force(r0[i], r0[j]);
                for(int k=0; k<ndim; k++){
                    r1[i][k] = r0[i][k] + v0[i][k]*dt + fij[k] * force_coefficient;
                    r1[j][k] = r0[j][k] + v0[j][k]*dt - fij[k] * force_coefficient;
                }
            } 
        }
    }

};

int main(){
    dynamics md(10, 0.01, 10);
    md.initialize_positions_velocities();
    md.write_positions_to_file();
    for(int i=0; i<md.npart; i++){
        for(int j=0; j<ndim; j++){
            md.r0[i][j] = md.r1[i][j];
        }
    }
    md.step += 1;
    md.write_positions_to_file();
}

