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

    particle p1;
    particle p2;
    particle p3;

    p1.pos = {1, 1, 1};
    p1.dir = {0, 0, 1};

    p2.pos = {2, 1, 1};
    p2.dir = {0, 0, 1};

    p3.pos = {1, 1, 2};
    p3.dir = {0, 1, 0};

    cout << "p1" << endl;
    print(p1);

    cout << "p2" << endl;
    print(p2);

    cout << "p3" << endl;
    print(p3);

    cout << endl;

    p2 = p3;

    p2.pos[0] += 1;

    cout << "p1" << endl;
    print(p1);

    cout << "p2" << endl;
    print(p2);

    cout << "p3" << endl;
    print(p3);
}
