#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>

#define ndim 3

using namespace std;

#include "tools.hpp"
#include "simulation.hpp"

void distance_test(){
    // test to see if minimum distance is working

    particle p1;
    particle p2;
    double box_l = 10;
    double L = 1;

    p1.pos = {1, 1, 1};
    p1.dir = {0, 0, 1};

    p2.pos = {2, 1, 1};
    p2.dir = {0, 0, 1};

    // this should be 1
    cout << dist_rods(p1, p2, box_l, L) << endl;

    p1.pos = {1, 1, 2};
    p1.dir = {0, 1, 0};

    p2.pos = {1, 1, 1};
    p2.dir = {0, 0, 1};

    // this should be 0.5
    cout << dist_rods(p1, p2, box_l, L) << endl;

    p1.pos = {1, 1, 1};
    p1.dir = {1, 0, 0};

    p2.pos = {3, 1, 1};
    p2.dir = {1, 0, 0};

    // this should be 0.5
    cout << dist_rods(p1, p2, box_l, L) << endl;
}

void write_to_file_test(){
    particle p1;
    particle p2;
    double box_l = 10;
    double L = 1;
    double D = 2;

    p1.pos = {1, 1, 2};
    p1.dir = {0, 1, 0};

    p2.pos = {1, 1, 1};
    p2.dir = {0, 0, 1};

    ofstream file;
    file.open("coords.dat");

    init_particle_file(file, 2, box_l);

    write_particle_to_file(file, p1, L, D);
    write_particle_to_file(file, p2, L, D);

    file.close();
}

int main(){
    write_to_file_test();
}
