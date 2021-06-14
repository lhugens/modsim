struct simulation{
    int N;                              // number of particles
    int step = 0;                       // current step of simulation
    int accept = 0;                     // number of accepted proposals
    double L;                           // length of the spherocylinders
    double D  = 1;                      // Diameter of the cylinders
    double D2 = pow(D, 2);              // Diameter**2
    double dl = 0.02;                   // maximum proposed displacement for each component
    double dn = 0.05;                   // maximum proposed change in direction for each component
    double dd = 0.001;                  // maximum proposed change in rho
    double rho;                         // density of the system (volume occupied by particles / total volume)
    double v0;                          // volume of one particle
    double rho_cp;                      // close packed density 
    double accept_rate;                 // acceptance rate
    vector<double> box_l;               // coordinates of box corners, one is at the origin
    vector<double> box_l_proposed;      // coordinates of box corners, one is at the origin
    vector<particle> part;              // main particle matrix
    vector<particle> part_proposed;     // main particle matrix

    simulation(double l) : L(l) 
    {
        v0 = M_PI * (L * pow(D, 2) / 4 + pow(D, 3) / 6);
        rho_cp = 2 / (sqrt(2) + (L/D)*sqrt(3));
    };

    inline void write_config(int step_no){write_config_to_file(step_no, N, box_l, part, L, D);}

    void fcc_config(int N_side, double rho_initial){
        auto [part_c, box_l_c, N_c] = fcc_config_initial(N_side, L);
        part  = part_c;
        box_l = box_l_c;
        N     = N_c;
        part_proposed = part;
        box_l_proposed = box_l;
        write_config(0);
        set_rho(rho_initial);
        part = part_proposed;
        box_l = box_l_proposed;
        //cout << "rho " << N * v0 / (box_l[0] * box_l[1] * box_l[2] * rho_cp) << endl;
    }

    void set_rho(double r){
        rho = r;

        double current_volume = box_l[0] * box_l[1] * box_l[2];
        double desired_volume = N * v0 / (rho * rho_cp);
        double scale_factor   = cbrt(desired_volume / current_volume);

        for(int j=0; j<ndim; j++){
            for(int i=0; i<N; i++){
                part_proposed[i].pos[j] *= scale_factor;
            }
            box_l_proposed[j] *= scale_factor;
        }
    }

    void propose_dl_dn(){
        step++;
        for(int i=0; i<N; i++){
            for(int j=0; j<ndim; j++){
                part_proposed[i].pos[j] = part[i].pos[j] + (2 * rdouble() - 1) * dl;
                part_proposed[i].dir[j] = part[i].dir[j] + (2 * rdouble() - 1) * dn;

                // boundary conditions
                part_proposed[i].pos[j] -= box_l[j] * floor(part_proposed[i].pos[j]/box_l[j]);
            }
            normalize(part_proposed[i].dir);
        }
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
        print(p1);
        print(p2);
        cout << sqrt(dist_rods(p1, p2, box_l, L)) << endl;
        cout << endl;
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
        if(!exists_overlap()){
            accept_proposal();
            accept++;
        }
        accept_rate = (double)accept/step;
    }

};


