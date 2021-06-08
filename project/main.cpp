#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>

#define ndim 3

using namespace std;

#include "mt19937.h"
#include "tools.hpp"
#include "simulation.hpp"

int main(){
    dsfmt_seed(time(NULL));

    simulation s;
    s.fcc_config(4, 0.2);
    cout << s.exists_overlap() << endl;
    s.write_config(0);


    for(int step=1; step<25; step++){
        s.propose();
        s.metropolis_acceptance();
        s.write_config(step);
    }
}
