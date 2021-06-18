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

/*
void test_file_loader(){
    double D = 1.0;
    string folder = "";
    auto [part_c, box_l_c, N_c, L_c] = file_config("liquid_config.dat");
    write_config_to_file(0, N_c, box_l_c, part_c, L_c, D, folder);
}
*/


/*
tuple<vector<particle>, vector<double>, int> liquid_config_initial(int N_part, double L, double rho){
    vector<particle> part;
    vector<double> box_l(ndim);
    particle p_temp;
    double side = 2*L;
    double temp = rdouble();

    box_l[0] = side*N_part + side;
    box_l[1] = side*N_part + side;
    box_l[2] = side*N_part + side;

    for(int i = 0; i < N_part + 1; i++){
        for(int j = 0; j < N_part + 1; j++){
            for(int l = 0; l < N_part + 1; l++){
                p_temp.pos[0] = (double)(i+rdouble()*0.65)*side;
                p_temp.pos[1] = (double)(j+rdouble()*0.65)*side;
                p_temp.pos[2] = (double)(l+rdouble()*0.65)*side;
                double val = rdouble();
                if(val < rdouble()){
                    //cout << "change x dir" << endl;
                    p_temp.dir[0] = temp*0.5;
                    p_temp.dir[1] = 0.5;     
                    p_temp.dir[2] = 0.5;
                }
                if(val < rdouble()){
                    //cout << "change y dir" << endl;
                    p_temp.dir[0] = 0.5;
                    p_temp.dir[1] = temp*0.5;
                    p_temp.dir[2] = 0.5; 
                }
                if(val < rdouble()){
                    //cout << "change z dir" << endl;
                    p_temp.dir[0] = 0.5;
                    p_temp.dir[1] = 0.5; 
                    p_temp.dir[2] = temp*0.5;
                } 
    	        part.push_back(p_temp);
            }
        }   
    }     

    double D = 1;
    int N = part.size();
    double v0 = M_PI * (L * pow(D, 2) / 4 + pow(D, 3) / 6);
    double rho_cp = 2 / (sqrt(2) + (L/D)*sqrt(3));
    double current_volume = box_l[0] * box_l[1] * box_l[2];
    double desired_volume = N * v0 / (rho * rho_cp);
    double scale_factor   = sqrt(desired_volume / current_volume);

    for(int j=0; j<ndim-1; j++){
        for(int i=0; i<N; i++){
            part[i].pos[j] *= scale_factor;
        }
        box_l[j] *= scale_factor;
    }
    return {part, box_l, N};
}

tuple<vector<particle>, vector<double>, int> liquid_config_initial(int N_part, double L, double rho){
    vector<particle> part;
    vector<double> box_l(ndim);
    double side = 2*L;
    double temp = rdouble();
    int N = 666;
    part.resize(N);

    box_l[0] = side*N_part + side;
    box_l[1] = side*N_part + side;
    box_l[2] = side*N_part + side;

    for(int i=0; i<N; i++){
        for(int j=0; j<ndim; j++){
            part[i].pos[j] = box_l[j] * rdouble();
            part[i].dir[j] = (2 * rdouble() - 1);
            normalize(part[i].dir);
        }
    }

    double D = 1;
    double v0 = M_PI * (L * pow(D, 2) / 4 + pow(D, 3) / 6);
    double rho_cp = 2 / (sqrt(2) + (L/D)*sqrt(3));
    double current_volume = box_l[0] * box_l[1] * box_l[2];
    double desired_volume = N * v0 / (rho * rho_cp);
    double scale_factor = cbrt(desired_volume / current_volume);

    for(int j=0; j<ndim; j++){
        for(int i=0; i<N; i++){
            part[i].pos[j] *= scale_factor;
        }
        box_l[j] *= scale_factor;
    }
    return {part, box_l, N};
}

void liquid_config(int N_side, double rho_initial){
    bool overlap = true;
    cout << "trying to find liquid phase" << endl;
    //while(overlap){
        rho = rho_initial;
        auto [part_c, box_l_c, N_c] = liquid_config_initial(N_side, L, rho_initial);
        part  = part_c;
        box_l = box_l_c;
        N     = N_c;
        V = box_l[0] * box_l[1] * box_l[2];
        overlap = exists_general_overlap();
    //}
    cout << "found liquid phase" << endl;
}

// run this to see that for our initial fcc, S = 1
void order_parameter(){
    simulation s("dont_matter", 3);
    s.fcc_config(5, 0.5);
    s.update_S();
}

NVT("coords_LIQUID/coords_", 0.5, 3.5, 1e5, 1e4);
NPT(0, "coords_EOS/coords_", 2.0, 0.3, 3, 1e3, 1, 0);
EOS();
NVT("LIQUID_CONFIG/coords", 0.3, 3.0, 1e5, 1e3);
*/

int main(){
    dsfmt_seed(time(NULL));

    vector<double> box_l = {10, 10, 10};
    vector<particle> part;

    particle p;

    p.pos = {1, 1, 1};
    p.dir = {0, 0, 1};

    part.push_back(p);

    p.pos = {3, 3, 3};
    p.dir = {0, 1, 0};

    part.push_back(p);

    p.pos = {5, 5, 5};
    p.dir = {1, 0, 0};

    part.push_back(p);

    string folder = "TEST/test"; 
    double L = 2;
    double D = 1;
    int N = 3;
    write_config_to_file(0, N, box_l, part, L, D, folder);
}
