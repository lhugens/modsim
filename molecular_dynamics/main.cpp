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
    double force_coefficient;   // force coefficient
    vector<vector<double>> r;  // array to store particle positions at t
    vector<vector<double>> dr;  // array to store particle displacement
    vector<vector<double>> f;   // array to store total force acting on each particle


    dynamics(int NPART, double DT, double BOX_L) :
        npart(NPART), dt(DT), box_l(BOX_L)
    {
        force_coefficient = pow(dt, 2) / m;
        r.resize(npart, vector<double>(ndim));
        dr.resize(npart, vector<double>(ndim));
        f.resize(npart, vector<double>(ndim));
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
                file << setw(20) << r[i][j];
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

    void update_forces(){
        double sq, factor;
        vector<double> rij(ndim);

        for(int i=0; i<npart-1; i++){
            for(int j=i+1; j<npart; j++){
                for(int k=0; k<ndim; k++){
                    rij[k]  = r[i][k]-r[j][k];
                    rij[k] -= round(rij[k] / box_l) * box_l; // periodic boundary conditions
                    f[i][k] = 0;
                    f[j][k] = 0;
                }
                sq = square_norm(rij);
                factor = 48 * (1/pow(sq, 3) - 0.5) / pow(sq, 4);
                for(int k=0; k<ndim; k++){
                    f[i][k] += rij[k] * factor;
                    f[j][k] -= rij[k] * factor;
                }
            }
        }
    }

    void initialize_positions_velocities(){

        // random r0 positions, each component in [-box_l, box_l]
        for(int i=0; i<npart; i++){
            for(int j=0; j<ndim; j++){
                r[i][j] = box_l * rand();
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
        update_forces();

        for(int i=0; i<npart; i++){
            for(int j=0; j<ndim; j++){
                dr[i][j] = v0[i][j]*dt + f[i][j] * force_coefficient / 2;
                r[i][j] += dr[i][j];
                r[i][j] += box_l - box_l * (int)((box_l + r[i][j])/box_l);
            }
        }
    }

    void verlet_step(){
        step += 1;
        update_forces();
        double temp;

        for(int k=0; k<npart; k++){
            for(int l=0; l<ndim; l++){
                dr[k][l] += f[k][l] * force_coefficient;
                r[k][l] += dr[k][l];
                r[k][l] += box_l - box_l * (int)((box_l + r[k][l])/box_l);
            }
        }
    }

};

int main(){
    dynamics md(10, 0.01, 10);
    md.initialize_positions_velocities();
    md.write_positions_to_file();
    for(int i=0; i<1000; i++){
        md.verlet_step();
        md.write_positions_to_file();
    }
}