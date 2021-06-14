struct simulation{
    int N;                              // number of particles
    int step = 0;                       // current step of simulation
    int accept = 0;                     // number of accepted proposals
    int pmoved;
    double L;                           // length of the spherocylinders
    double D  = 1;                      // Diameter of the cylinders
    double D2 = pow(D, 2);              // Diameter**2
    double dl = 0.02;                   // maximum proposed displacement for each component
    double dn = 0.05;                   // maximum proposed change in direction for each component
    double drho = 0.001;                // maximum proposed change in rho
    double betaP = 1;
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

    // random integers from in the interval [0, N-1] 
    inline double rint(){return (int)(dsfmt_genrand() * N);}


    void fcc_config(int N_side, double rho_initial){
        auto [part_c, box_l_c, N_c] = fcc_config_initial(N_side, L, rho_initial);
        part  = part_c;
        box_l = box_l_c;
        N     = N_c;
        part_proposed = part;
        box_l_proposed = box_l;
        write_config(0);
        cout << "rho " << N * v0 / (box_l[0] * box_l[1] * box_l[2] * rho_cp) << endl;
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

    void propose_NVT(){
        pmoved = rint();
        if(rand() < 0.5){
            for(int j=0; j<ndim; j++){
                part_proposed[pmoved].pos[j] = part[pmoved].pos[j] + (2 * rdouble() - 1) * dl;
                part_proposed[pmoved].pos[j] -= box_l[j] * floor(part_proposed[pmoved].pos[j]/box_l[j]);
            }
        }
        else{
            for(int j=0; j<ndim; j++){
                part_proposed[pmoved].dir[j] = part[pmoved].dir[j] + (2 * rdouble() - 1) * dn;
            }
            normalize(part_proposed[pmoved].dir);
        }
    }

    void propose_NPT(){
        propose_NVT();
        set_rho(rho + drho * (2 * rdouble() - 1));
    }

    void accept_proposal(){
        accept++;
        accept_rate = (double)accept/step;
        for(int i=0; i<N; i++){
            part[i] = part_proposed[i];
        }
        for(int i =0; i<ndim; i++){
            box_l[i] = box_l_proposed[i];
        }
    }

    bool they_overlap(particle &p1, particle &p2){
        //print(p1);
        //print(p2);
        //cout << sqrt(dist_rods(p1, p2, box_l, L)) << endl;
        //cout << endl;
        return (dist_rods(p1, p2, box_l, L) < D2);
    }

    bool exists_initial_overlap(){
        for(int i=0; i<N; i++){
            for(int j=1; j<N-1; j++){
                if(they_overlap(part_proposed[i], part_proposed[j]))
                    return true;
            }
        }
        return false;
    }

    bool exists_overlap(){
        for(int i=0; i<N; i++){
            if(pmoved != i){
                if(they_overlap(part_proposed[pmoved], part_proposed[i]))
                    return true;
            }
        }
        return false;
    }

    void metropolis_acceptance_NVT(){
        if(!exists_overlap()){
            accept_proposal();
        }
    }

    void metropolis_acceptance_NPT(){
        if(!exists_overlap()){
            double V      = pow(box_l[0], 3);
            double V_prop = pow(box_l_proposed[0], 3);
            if(rdouble() < exp(-betaP*(V_prop - V) + N * log(V_prop/V))){
                accept_proposal();
            }
        }
    }

};


