#include <iostream>
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
    double beyta = 0.5;                   // beyta = 1 / (k_B * T)
    double nu;                            // frequency of stochastic collisions
    double andersen_probability;          // frequency of stochastic collisions
    double force_coefficient;             // force coefficient
    double kinetic_energy;                // total kinetic_energy of the system
    double potential_energy;              // total potential_energy of the system
    int average_collisions_bath = 0.0;    // number of average collisions woth heat bath
    vector<vector<double>> r;             // array to store particle positions at t
    vector<vector<double>> dr;            // array to store particle displacement
    vector<vector<double>> f;             // array to store total force acting on each particle
    default_random_engine gen;            // random_seed 
    normal_distribution<double> gaussian; // gaussian distribution
    ofstream file;

    dynamics(int NPART, double DT, double BOX_L, double NU) :
        npart(NPART), dt(DT), box_l(BOX_L), nu(NU)
    {
        force_coefficient = pow(dt, 2) / m;
        r.resize(npart, vector<double>(ndim));
        dr.resize(npart, vector<double>(ndim));
        f.resize(npart, vector<double>(ndim));
        dsfmt_seed(time(NULL));
        file.open("state_variables.txt");

        r_cut_squared = pow(box_l / 3, 2);
        double r_cut_pow_6 = pow(r_cut_squared, 3);
        e_cut = (4 / r_cut_pow_6) * ((1/r_cut_pow_6) - 1);
        double standard_deviation = 1 / sqrt(m*beyta);
        gaussian = normal_distribution<double>(0.0, standard_deviation);
        andersen_probability = nu*dt;
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
            file << setw(20) << 0.5 << endl;
        }
        file.close();
    }

    void write_state_variables_to_file(){
        update_kinetic_energy();
        file << left << setw(20) << step 
                     << setw(20) << kinetic_energy 
                     << setw(20) << potential_energy
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
        double sq, factor;
        potential_energy = 0;
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
                if (sq < r_cut_squared){
                    factor = 48 * (1/pow(sq, 3) - 0.5) / pow(sq, 4);
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

        write_positions_to_file();

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
        update_forces();

        step += 1;

        for(int i=0; i<npart; i++){
            for(int j=0; j<ndim; j++){
                dr[i][j] = v0[i][j]*dt + f[i][j] * force_coefficient / 2;
                r[i][j] += dr[i][j];
                r[i][j] += box_l - box_l * (int)((box_l + r[i][j])/box_l);
            }
        }
    }

    void update_kinetic_energy(){
        kinetic_energy = 0;
        for(int k=0; k<npart; k++){
            for(int l=0; l<ndim; l++){
                kinetic_energy += m * pow((dr[k][l] / dt), 2) / 2;
            }
        }
    }

    void verlet_step_NVE(){
        step += 1;
        update_forces();

        for(int k=0; k<npart; k++){
            for(int l=0; l<ndim; l++){
                dr[k][l] += f[k][l] * force_coefficient;
                r[k][l] += dr[k][l];
                r[k][l] += box_l - box_l * (int)((box_l + r[k][l])/box_l);
            }
        }
    }

    void verlet_step_NVT(){
        step += 1;
        update_forces();

        for(int k=0; k<npart; k++){
            if(rand() < andersen_probability){
                average_collisions_bath += 1;
                for(int l=0; l<ndim; l++){
                    dr[k][l] = boltzmann_velocity() * dt;
                    r[k][l] += dr[k][l];
                    r[k][l] += box_l - box_l * (int)((box_l + r[k][l])/box_l);
                }
            }
            else{
                for(int l=0; l<ndim; l++){
                    dr[k][l] += f[k][l] * force_coefficient;
                    r[k][l] += dr[k][l];
                    r[k][l] += box_l - box_l * (int)((box_l + r[k][l])/box_l);
                }
            }
        }
    }

};

void test_verlet(){
    dynamics md(2, 0.00001, 5, 1);

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
    md.dr[0][0] =  0.01;
    md.dr[1][0] = -0.01;
    for(int i=0; i<md.npart; i++){
        for(int j = 1; j<ndim; j++){
            md.dr[i][j] = 0;
        }
    }

    cout << "particle dr's:" << endl;
    md.print(md.dr);
    md.update_forces();
    md.print(md.f);

    int total_steps = 1000;
    md.write_positions_to_file();
    for(int i=0; i<total_steps; i++){
        md.verlet_step_NVE();
        md.write_positions_to_file();
        md.write_state_variables_to_file();
        cout << "\r [" << setw(3) << round((double)i * 100 /total_steps)  << "%] " << flush;
    }
}

void run_NVE(){
    dynamics md(100, 0.00001, 5, 1);
    md.initialize_positions_velocities();
    md.write_positions_to_file();
    int total_steps = 1000;
    for(int i=2; i<total_steps; i++){
        md.verlet_step_NVE();
        md.write_positions_to_file();
        md.write_state_variables_to_file();
        cout << "\r [" << setw(3) << round((double)i * 100 /total_steps)  << "%] " << flush;
    }
}

void run_NVT(){
    dynamics md(100, 0.0001, 5, 10);
    md.initialize_positions_velocities();
    md.write_positions_to_file();
    int total_steps = 1000;
    for(int i=2; i<total_steps; i++){
        md.verlet_step_NVT();
        md.write_positions_to_file();
        md.write_state_variables_to_file();
        cout << "\r [" << setw(3) << round((double)i * 100 /total_steps)  << "%] " << flush;
    }
}

int main(){
    //test_verlet();
    //run_NVE();
    run_NVT();
}










