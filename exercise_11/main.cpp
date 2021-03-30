#include <stdio.h>
#include <vector>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 512

// Initialization variables
const int    mc_steps      = 10000;
const int    output_steps  = 100;
double       density       = 0.8;
double       mg            = 5.5;
double       dz            = 0.1;
double       beta          = 0.5;
const double delta         = 0.1;
const double r_cut         = 2.5;
//const char*  init_filename = "fcc.dat";

// Simulation variables
int n_particles = N;    // number of particles
double e_cut;           // cut energy
double box_l = 8.0; 
double r[N][NDIM];      // coordinates j of particle i
double box[NDIM] = {box_l, box_l, box_l};
double vol;             // volume of the box

double energy = 0.0;    
double virial = 0.0;

//function generates a random initial configuration
//of N particles inside box

void generate_data(){
    for(int p=0; p<N; p++){
        for(int i=0; i<NDIM; i++){
            r[p][i] = dsfmt_genrand() * box[i];
        }
    }
}

// function writes the current positions of the particles 
// to a file, in the format:
// 500
// 0.000000 8.549880
// 0.000000 8.549880
// 0.000000 8.549880
// 0.865781        0.913697        0.792968        1.000000
// 1.709977        1.709977        0.854988        1.000000
// 0.854988        1.709977        1.709977        1.000000

void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coords/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", 1.0);
    }
    fclose(fp);
}

//function creates a histogram of the z values
//of all the all the particles

void write_histogram(int index, int step){
    char buffer[128];
    sprintf(buffer, "hist/hist_mg%03d_step%07d.dat", index, step);
    FILE* fp = fopen(buffer, "w");
    std::vector<int> histogram;
    int bin;

    for(int p=0; p<n_particles; p++){
        bin = int(r[p][NDIM-1]/dz);
        if (bin > static_cast<int>(histogram.size()) - 1)
            histogram.resize(bin + 1);

        histogram[bin] += 1;
    }

    int i = 0;
    for(auto & value : histogram){
        fprintf(fp, "%lf %d\n",(double)i*dz,value);
        i++;
    }
    fclose(fp);
}

// particle_info_t struct
typedef struct{
    double energy = 0.0;
    double virial = 0.0;
}particle_info_t;

// measurement_t struct
typedef struct{
    double average_pressure = 0.0;
    double mu_excess = 0.0;
}measurement_t;

// particle_energy_and_virial(pid) alculates 
// the energy contribution 
// and the contrib. to the virial term of a single 
// particle chosen at random

particle_info_t particle_energy_and_virial(int pid){

    particle_info_t info;
    for(int n = 0; n < n_particles; n++){
        if(n == pid) 
            continue;

        double dist2 = 0.0;
        for(int d = 0; d < NDIM; d++){
            double min_d = r[pid][d] - r[n][d];
            min_d -= (int)(2.0 * min_d / box[d]) * box[d];
            dist2 += pow(min_d, 2);
        }

        if(dist2 <= pow(r_cut, 2)){
            double temp = 1.0 / pow(dist2, 3);
            info.energy += 4.0 * temp * (temp - 1.0) - e_cut;
            info.virial += 24.0 * temp * (2.0 * temp - 1.0);
        }
    }
    info.energy += mg*r[pid][NDIM-1];
    return info;
}

// measure() calculates the average pressure
// by iterating through all particles and calling
// particle_energy_and_virial for each, and averaging over

measurement_t measure(){
    measurement_t result;
    result.average_pressure = density/beta * virial / (3 * vol);
    return result;
}

// move_particle() attempts to displace a single random
// particle, compares it's old and new contributions to the total energy
// and to the virial term, and accepts/rejects (returns 1/0) according
// to Metropolis

int move_particle(){
    int rpid = n_particles * dsfmt_genrand();

    particle_info_t info = particle_energy_and_virial(rpid);

    double old_pos[NDIM];

    // x and y directions behave normally
    for(int d = 0; d < NDIM-1; d++){
        old_pos[d] = r[rpid][d];
        r[rpid][d] += delta * (2.0 * dsfmt_genrand() - 1.0) + box[d];
        r[rpid][d] -= (int)(r[rpid][d] / box[d]) * box[d];
    }

    // z direction
    old_pos[NDIM-1] = r[rpid][NDIM-1];
    r[rpid][NDIM-1] += delta * (2.0 * dsfmt_genrand() - 1.0);

    if(r[rpid][NDIM-1] <= 0){
        for(int d = 0; d < NDIM; d++) 
            r[rpid][d] = old_pos[d];
        return 0;
    }

    particle_info_t new_info = particle_energy_and_virial(rpid);

    double dE = new_info.energy - info.energy;
    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        energy += dE;
        virial += new_info.virial - info.virial;
        return 1;
    }

    for(int d = 0; d < NDIM; d++) 
        r[rpid][d] = old_pos[d];

    return 0;
}



// set_density() takes any array of position r
// and correspoding box, and scales both in order
// to achieve a desired density

void set_density(){
    double volume = 1.0;
    for(int d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = n_particles / density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(int n = 0; n < n_particles; ++n){
        for(int d = 0; d < NDIM; ++d) 
            r[n][d] *= scale_factor;
    }
    vol = 1.0;
    for(int d = 0; d < NDIM; ++d){
        box[d] *= scale_factor;
        vol *= box[d];
    }
}

int main(int argc, char* argv[]){

    assert(delta > 0.0);

    e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0));


    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    size_t seed = time(NULL);
    dsfmt_seed(seed);

    double mgs[] = {1, 2.5, 5};

    // Main loop

    for(int mg_i = 0; mg_i < 3; mg_i++){
        mg = mgs[mg_i];

        // Initial configuration
        
        generate_data();
        
        for(int d = 0; d < NDIM; ++d) 
            assert(r_cut <= 0.5 * box[d]);

        for(int n = 0; n < n_particles; ++n){
            particle_info_t info = particle_energy_and_virial(n);
            energy += info.energy;
            virial += info.virial;
        }
        energy *= 0.5;
        virial *= 0.5;

        double volume = 1.0;
        for(int d = 0; d < NDIM; ++d) 
            volume *= box[d];

        //printf("Starting volume: %f\n",  volume);
        //printf("Starting energy: %f\n",  energy);
        //printf("Starting virial: %f\n",  virial);
        //printf("Starting seed:   %lu\n", seed);
        printf("Density:         %f\n",  density);
        printf("Beta:            %f\n",  beta);

        // Main procedure
        
        int accepted = 0;
        for(int step = 0; step < mc_steps; ++step){
            for(int n = 0; n < n_particles; ++n){
                accepted += move_particle();
            }

            if(step % output_steps == 0){
                printf("Step %d. Move acceptance: %f.\n",
                    step, (double)accepted / (n_particles * output_steps)
                );
                accepted = 0;
                write_data(step);
                write_histogram(mg_i, step);
            }
        }
    }

    return 0;
}
