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
    s.fcc_config(4, 0.3);
    bool overlap = s.exists_overlap();
    cout << overlap << endl;
    assert(!overlap);
    s.write_config(0);


    for(int step=1; step<1000; step++){
        s.propose();
        s.metropolis_acceptance();
        cout << s.density << endl;
        s.write_config(step);
    }
}
