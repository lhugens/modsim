#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"
//#pragma warning(disable : 5208)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 512


// Initialization variables
const int    mc_steps = 500;
const int    output_steps = 100;
double       density = 0.8;
double       beta = 0.5;
const double delta = 0.1;
const double r_cut = 2.5;
const char* init_filename = "fcc.dat";

// Simulation variables
int n_particles = 0;    // number of particles
double e_cut;           // cut energy
double r[N][NDIM];      // coordinates j of particle i
double box[NDIM];       // box max_i
double vol;          // volume of the box

double energy = 0;
double virial = 0;

// function reads the initial particle configuration 
// from an external file. Sets the global variables
// r and box 

void read_data() {
    FILE* fp = fopen(init_filename, "r");
    double dmin, dmax;
    fscanf(fp, "%d\n", &n_particles);
    for (int d = 0; d < NDIM; ++d) {
        fscanf(fp, "%lf %lf\n", &dmin, &dmax);
        box[d] = fabs(dmax - dmin);
    }
    for (int n = 0; n < n_particles; ++n) {
        for (int d = 0; d < NDIM; ++d)
            fscanf(fp, "%lf\t", &r[n][d]);
        double diameter;
        fscanf(fp, "%lf\n", &diameter);
    }
    fclose(fp);
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

void write_data(int step) {
    char buffer[600];
    sprintf(buffer, "coords/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for (d = 0; d < NDIM; ++d) {
        fprintf(fp, "%lf %lf\n", 0.0, box[d]);
    }
    for (n = 0; n < n_particles; ++n) {
        for (d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", 1.0);
    }
    fclose(fp);
}

// particle_info_t struct
typedef struct {
    double energy;
    double virial;
}particle_info_t;

// measurement_t struct
typedef struct {
    double average_pressure;
    double mu_excess;
}measurement_t;


// particle_energy_and_virial(pid) alculates 
// the energy contribution 
// and the contrib. to the virial term of a single 
// particle chosen at random

particle_info_t particle_energy_and_virial(int pid) {

    particle_info_t info;
    info.energy = 0.0;
    info.virial = 0.0;

    for (int n = 0; n < n_particles; n++) {
        if (n == pid)
            continue;

        double dist2 = 0.0;
        for (int d = 0; d < NDIM; d++) {
            double min_d = r[pid][d] - r[n][d];
            min_d -= (int)(2.0 * min_d / box[d]) * box[d];
            dist2 += pow(min_d, 2);
        }

        if (dist2 <= pow(r_cut, 2)) {
            double temp = 1.0 / pow(dist2, 3);
            info.energy += 4.0 * temp * (temp - 1.0) - e_cut;
            info.virial += 24.0 * temp * (2.0 * temp - 1.0);
        }
    }
    return info;
}

// measure() calculates the average pressure
// by iterating through all particles and calling
// particle_energy_and_virial for each, and averaging over6

measurement_t measure() {
    measurement_t result;
    result.average_pressure = density / beta + virial / (3 * vol);// leo: density / beta * virial / (3 * vol);
    result.mu_excess =  log(density)/beta + energy/(beta*n_particles);//mu ideal + mu excess
    return result;
}//ideal gas = k T ln(rho * lambda)

// move_particle() attempts to displace a single random
// particle, compares it's old and new contributions to the total energy
// and to the virial term, and accepts/rejects (returns 1/0) according
// to Metropolis

int move_particle() {
    int rpid = n_particles * dsfmt_genrand();

    particle_info_t info = particle_energy_and_virial(rpid);

    double old_pos[NDIM];
    for (int d = 0; d < NDIM; d++) {
        old_pos[d] = r[rpid][d];
        r[rpid][d] += delta * (2.0 * dsfmt_genrand() - 1.0) + box[d];
        r[rpid][d] -= (int)(r[rpid][d] / box[d]) * box[d];
    }

    particle_info_t new_info = particle_energy_and_virial(rpid);

    double dE = new_info.energy - info.energy;
    if (dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)) {
        energy += dE;
        virial += new_info.virial - info.virial;
        return 1;
    }

    for (int d = 0; d < NDIM; d++)
        r[rpid][d] = old_pos[d];

    return 0;
}

// set_density() takes any array of position r
// and correspoding box, and scales both in order
// to achieve a desired density

void set_density() {
    double volume = 1.0;
    for (int d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = n_particles / density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for (int n = 0; n < n_particles; ++n) {
        for (int d = 0; d < NDIM; ++d)
            r[n][d] *= scale_factor;
    }
    vol = 1.0;
    for (int d = 0; d < NDIM; ++d) {
        box[d] *= scale_factor;
        vol *= box[d];
    }
}

int main(int argc, char* argv[]) {

    assert(delta > 0.0);

    e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0));

    read_data();

    if (n_particles == 0) {
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    size_t seed = time(NULL);
    dsfmt_seed(seed);

    // Build density and temperature array

    const int len_rhos = 10;
    double min_rho = 0.1;
    double max_rho = 1.1;
    double rhos[len_rhos];
    double d_rho = (max_rho - min_rho) / (double)len_rhos;

    double betas[] = { 0.5, 1.0, 2.0 };


    // Main loop

    for (int beta_i = 0; beta_i < 3; beta_i++) {
        beta = betas[beta_i];
        for (int rho_i = 0; rho_i < len_rhos; rho_i++) {
            double temp = rho_i;
            density = (temp + 1) * d_rho;
            rhos[rho_i] = density;

            // Initial configuration

            set_density();

            for (int d = 0; d < NDIM; ++d)
                assert(r_cut <= 0.5 * box[d]);

            for (int n = 0; n < n_particles; ++n) {
                particle_info_t info = particle_energy_and_virial(n);
                energy += info.energy;
                virial += info.virial;
            }
            energy *= 0.5;
            virial *= 0.5;

            double volume = 1.0;
            for (int d = 0; d < NDIM; ++d)
                volume *= box[d];

            //printf("Starting volume: %f\n",  volume);
            //printf("Starting energy: %f\n",  energy);
            //printf("Starting virial: %f\n",  virial);
            //printf("Starting seed:   %lu\n", seed);
            printf("Density:         %f\n", density);
            printf("Beta:            %f\n", beta);
            char buffer[600];
            sprintf(buffer, "measurements_beta_%03d_rho_%03d.dat", beta_i, rho_i);

            FILE* fp = fopen(buffer, "w");

            // Main procedure

            int accepted = 0;
            for (int step = 0; step < mc_steps; ++step) {
                for (int n = 0; n < n_particles; ++n) {
                    accepted += move_particle();
                }

                measurement_t ms = measure();

                fprintf(fp, "%d\t%f\t%f\n", step, ms.average_pressure, ms.mu_excess);


                if (step % output_steps == 0) {
                    printf("Step %d. Move acceptance: %f.\n",
                        step, (double)accepted / (n_particles * output_steps)
                    );
                    accepted = 0;
                    //write_data(step);
                }
            }

            fclose(fp);

        }
    }

    return 0;
}