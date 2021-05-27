#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "tools.hpp"
#include "structs.hpp"

int main(){
    simulation s(10);

    // test to see if minimum distance is working

    particle p1;
    particle p2;

    /*
    p1.pos = {1, 1, 1};
    p1.dir = {0, 0, 1};

    p2.pos = {2, 1, 1};
    p2.dir = {0, 0, 1};

    // this should be 1
    cout << s.dist_rods(p1, p2) << endl;

    */

    p1.pos = {1, 1, 2};
    p1.dir = {0, 1, 0};

    p2.pos = {1, 1, 1};
    p2.dir = {0, 0, 1};

    // this should be 0.5
    cout << s.dist_rods(p1, p2) << endl;

    /*
    p1.pos = {1, 1, 1};
    p1.dir = {1, 0, 0};

    p2.pos = {3, 1, 1};
    p2.dir = {1, 0, 0};

    // this should be 0.5
    cout << s.dist_rods(p1, p2) << endl;
    */
}
