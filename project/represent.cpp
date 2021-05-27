#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

// trying to draw some particles

void write_particle(ofstream & file, double x, double y, double z, double m11, double m12, double m13, 
                                                                   double m21, double m22, double m23, 
                                                                   double m31, double m32, double m33){

    // particle type
    file << left << setw(20) << "SPHEROCYLINDER(0.5)";

    // particle center position
    file << setw(10) << x
         << setw(10) << y
         << setw(10) << z;

    // orientation
    file << setw(10) << m11 << setw(10) << m12 << setw(10) << m13;
    file << setw(10) << m21 << setw(10) << m22 << setw(10) << m23;
    file << setw(10) << m31 << setw(10) << m32 << setw(10) << m33;

    file << endl;
}

int main(){
    ofstream f;
    f.open("test.dat");

    // number of structures
    f << 2 << endl;

    // box dimensions
    f << left << setw(10) << 0.0 << setw(10) << 10.0 << endl;
    f << left << setw(10) << 0.0 << setw(10) << 10.0 << endl;
    f << left << setw(10) << 0.0 << setw(10) << 10.0 << endl;

    write_particle(f, 1, 1, 1, 1, 0, 0,
                               0, 1, 0,
                               0, 0, 1);

    write_particle(f, 3, 3, 3, 2, 0, 0,
                               0, 0.5, 0,
                               0, 0, 1);

    f.close();
    
}
