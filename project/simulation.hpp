struct simulation{
    int N;                              // number of particles
    double L  = 1;                      // length of the spherocylinders
    double D  = 1;                      // Diameter of the cylinders
    double D2 = pow(D, 2);              // Diameter**2
    double dl = 0.01;                    // maximum proposed displacement for each component
    double dn = 0.05;                     // maximum proposed change in direction for each component
    double dd = 0.001;                    // maximum proposed change in density
    double density;                     // density of the system (volume occupied by particles / total volume)
    vector<double> box_l;               // coordinates of box corners, one is at the origin
    vector<double> box_l_proposed;      // coordinates of box corners, one is at the origin
    vector<particle> part;              // main particle matrix
    vector<particle> part_proposed;     // main particle matrix

    inline void write_config(int step){write_config_to_file(step, N, box_l, part, L, D);}

    void fcc_config(int N_side, double density_initial){
        auto [part_c, box_l_c, N_c] = fcc_config_initial(N_side, L);
        part  = part_c;
        box_l = box_l_c;
        N     = N_c;
        part_proposed = part;
        box_l_proposed = box_l;
        set_density(density_initial);
        part = part_proposed;
        box_l = box_l_proposed;
    }

    void set_density(double dens){
        density = dens;

        double v0 = M_PI * (L * pow(D, 2) / 4 + pow(D, 3) / 6);
        double current_volume = box_l[0] * box_l[1] * box_l[2];
        double desired_volume = N * v0 / (2 * density/(sqrt(2) + (L/D)*sqrt(3)));
        double scale_factor   = cbrt(desired_volume / current_volume);

        for(int j=0; j<ndim; j++){
            for(int i=0; i<N; i++){
                part_proposed[i].pos[j] *= scale_factor;
            }
            box_l_proposed[j] *= scale_factor;
        }
    }

    void propose(){
        for(int i=0; i<N; i++){
            for(int j=0; j<ndim; j++){
                part_proposed[i].pos[j] = part[i].pos[j] + (2 * rdouble() - 1) * dl;
                part_proposed[i].dir[j] = part[i].dir[j] + (2 * rdouble() - 1) * dn;

                // boundary conditions
                part_proposed[i].pos[j] -= box_l[j] * floor(part_proposed[i].pos[j]/box_l[j]);
            }
            normalize(part_proposed[i].dir);
        }
        set_density(density + (2 * rdouble() - 1) * dd);
    }

    void accept_proposal(){
        for(int i=0; i<N; i++){
            part[i] = part_proposed[i];
        }
        for(int i =0; i<ndim; i++){
            box_l[i] = box_l_proposed[i];
        }
    }

    bool they_overlap(particle &p1, particle &p2){
        return (dist_rods(p1, p2, box_l, L) < D2);
    }

    bool exists_overlap(){
        for(int i=0; i<N-1; i++){
            for(int j=i+1; j<N; j++){
                if(they_overlap(part_proposed[i], part_proposed[j])){
                    return true;
                }
            }
        }
        return false;
    }

    void metropolis_acceptance(){
        if(!exists_overlap())
            accept_proposal();
    }

};


