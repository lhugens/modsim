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

    simulation s(10);
    s.fcc_config();
}
