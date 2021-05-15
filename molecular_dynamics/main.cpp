#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <random>

#include "mt19937.h"

using namespace std;

#define ndim 3

struct dynamics{
    int npart;                            // number of particles
    int step = 0;                         // current time step of the simulation
    double box_l;                         // (cubic) box side length
    double r_cut_squared;                 // cut-off square distance of the potential
    double e_cut;                         // cut-off energy of the potential
    double dt;                            // time step
    double m = 1;                         // particle mass
    double beyta = 2;                     // beyta = 1 / (k_B * T)
    double nu;                            // frequency of stochastic collisions
    double andersen_probability;          // probability of colliding with heat bath
    double diameter;                      // particle diameter
    double force_coefficient;             // force coefficient
    double kinetic_energy;                // total kinetic_energy of the system
    double potential_energy;              // total potential_energy of the system
    int average_collisions_bath = 0;      // number of average collisions woth heat bath
    vector<vector<double>> r;             // array to store particle positions at t
    vector<vector<double>> v;             // array to store particle displacement
    vector<vector<double>> f;             // array to store total force acting on each particle
    default_random_engine gen;            // random_seed 
    normal_distribution<double> gaussian; // gaussian distribution
    ofstream file;

    dynamics(int NPART, double DT, double BOX_L, double NU) :
        npart(NPART), dt(DT), box_l(BOX_L), nu(NU)
    {
        dsfmt_seed(time(NULL));

        r.resize(npart, vector<double>(ndim));
        v.resize(npart, vector<double>(ndim));
        f.resize(npart, vector<double>(ndim));

        file.open("state_variables.txt");

        force_coefficient = pow(dt, 2) / (2*m);

        andersen_probability = nu*dt;

        diameter = box_l / 20;

        r_cut_squared = pow(box_l / 3, 2);
        double r_cut_pow_6 = pow(r_cut_squared, 3);
        e_cut = (4 / r_cut_pow_6) * ((1/r_cut_pow_6) - 1);

        double standard_deviation = 1 / sqrt(m*beyta);
        gaussian = normal_distribution<double>(0.0, standard_deviation);

    }

    inline double rand(){
        return dsfmt_genrand();
    }

    inline double boltzmann_velocity(){
        return gaussian(gen);
    }

    void test_boltzmann_distribution(){
        ofstream file;
        file.open("boltzmann.txt");
        file << 1 / sqrt(m*beyta) << endl;

        for(int i=0; i<1000000; i++){
            file << boltzmann_velocity() << endl;
        }

        file.close();
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
            file << setw(20) << box_l/20 << endl;
        }
        file.close();
    }

    void write_state_variables_to_file(){
        update_kinetic_energy();
        file << left << setw(20) << step 
                     << setw(20) << kinetic_energy 
                     << setw(20) << potential_energy
                     << setw(20) << kinetic_energy + potential_energy
                     << setw(20) << (double)average_collisions_bath / (double)step
                     << endl;
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

    void update_forces(){
        potential_energy = 0;
        double sq, factor;
        vector<double> rij(ndim);

        for(int i=0; i<npart; i++){
            for(int j=0; j<ndim; j++){
                f[i][j] = 0.0;
            }
        }

        for(int i=0; i<npart-1; i++){
            for(int j=i+1; j<npart; j++){
                for(int k=0; k<ndim; k++){
                    rij[k]  = r[i][k]-r[j][k];
                    rij[k] -= round(rij[k] / box_l) * box_l; // periodic boundary conditions
                }

                sq = square_norm(rij);

                if (sq < r_cut_squared){
                    factor = 48 * ((1/pow(sq, 3)) - 0.5) / pow(sq, 4);
                    for(int k=0; k<ndim; k++){
                        f[i][k] += rij[k] * factor;
                        f[j][k] -= rij[k] * factor;
                    }
                    sq = pow(sq, 3); // now sq = r^6
                    potential_energy += (4 / sq) * ((1/sq) - 1);
                    potential_energy -= e_cut;
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

        update_forces();

        // randomize velocities, each component in [-1,1]
        // and compute total momentum
        vector<double> P(ndim);

        for(int j=0; j<ndim; j++){
            for(int i=0; i<npart; i++){
                v[i][j] = rand() * 2 - 1;
                P[j] += v[i][j];
            }
        }

        // shift velocities s.t. total momentum is 0
        for(int j=0; j<ndim; j++){
            for(int i=0; i<npart; i++){
                v[i][j] -= P[j] / (double)npart;
            }
        }

        // scale velocities s.t. E_kinetic = 3 N k_B T / 2
        double sum_square_speeds = 0;
        for(int i=0; i<npart; i++){
            sum_square_speeds += square_norm(v[i]);
        }
        double scaling_factor = sqrt(3 * npart / (beyta * m * sum_square_speeds));

        for(int i=0; i<npart; i++){
            for(int j=0; j<ndim; j++){
                v[i][j] *= scaling_factor;
            }
        }
    }
  
    void update_kinetic_energy(){
        kinetic_energy = 0;
        for(int k=0; k<npart; k++){
            for(int l=0; l<ndim; l++){
                kinetic_energy += m * pow((v[k][l]), 2) / 2;
            }
        }
    }

    void verlet_step_NVE(){
        step += 1;

        for(int k=0; k<npart; k++){
            for(int l=0; l<ndim; l++){
                r[k][l] += v[k][l] * dt + f[k][l] * force_coefficient;
                v[k][l] += f[k][l] * dt / (2*m);
                //r[k][l] += box_l - box_l * (int)((box_l + r[k][l])/box_l);
                r[k][l] -= box_l * floor(r[k][l]/box_l);
            }
        }
        update_forces();
        for(int k=0; k<npart; k++){
            for(int l=0; l<ndim; l++){
                v[k][l] += f[k][l] * dt / (2*m);
            }
        }
    }

    void thermostat(){
        for(int i=0; i<npart; i++){
            if(rand() < andersen_probability){
                average_collisions_bath++;
                for(int l=0; l<ndim; l++){
                    v[i][l] = boltzmann_velocity();
                }
            }
        }
    }
};

void test_verlet(){
    dynamics md(2, 0.001, 5, 1);

    // put the two particles on opposite ends of the box
    md.r[0][0] = 0;
    md.r[1][0] = md.box_l / 2;
    for(int i=0; i<md.npart; i++){
        for(int j = 1; j<ndim; j++){
            md.r[i][j] = 0;
        }
    }

    cout << "particle positions:" << endl;
    md.print(md.r);
    
    // particles going head on towards eachother
    md.v[0][0] =  1;
    md.v[1][0] = -1;
    for(int i=0; i<md.npart; i++){
        for(int j = 1; j<ndim; j++){
            md.v[i][j] = 0;
        }
    }

    cout << "particle v's:" << endl;
    md.print(md.v);
    md.update_forces();
    md.print(md.f);

    int total_steps = 10000;
    md.write_positions_to_file();
    for(int i=0; i<total_steps; i++){
        md.verlet_step_NVE();
        md.write_positions_to_file();
        md.write_state_variables_to_file();
        cout << "\r [" << setw(3) << round((double)i * 100 /total_steps)  << "%] " << flush;
    }
}

void run_NVE(int total_steps, int NPART, double DT, double BOX_L){
    dynamics md(NPART, DT, BOX_L, 1);
    md.initialize_positions_velocities();
    md.write_positions_to_file();
    md.verlet_step_NVE();
    for(int i=2; i<total_steps; i++){
        md.verlet_step_NVE();
        //md.write_positions_to_file();
        md.write_state_variables_to_file();
        cout << "\r [" << setw(3) << round((double)i * 100 /total_steps)  << "%] " << flush;
    }
}

void run_NVT(int total_steps, int NPART, double DT, double BOX_L, double NU){
    dynamics md(NPART, DT, BOX_L, NU);
    md.initialize_positions_velocities();
    md.write_positions_to_file();
    for(int i=1; i<total_steps; i++){
        md.verlet_step_NVE();
        md.thermostat();
        //md.write_positions_to_file();
        md.write_state_variables_to_file();
        cout << "\r [" << setw(3) << round((double)i * 100 /total_steps)  << "%] " << flush;
    }
}

int main(){
    //test_verlet();
    run_NVE(100000, 100, 0.001, 50);
    //run_NVT(100000, 100, 0.001, 50, 50);
}
