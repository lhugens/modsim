struct simulation{
    int N;                              // number of particles
    int step = 0;                       // current step of simulation
    int accept = 0;                     // number of accepted proposals
    int pmoved;                         // suggested particle to be displaced
    double L;                           // length of the spherocylinders
    double D  = 1;                      // diameter of the cylinders
    double D2 = pow(D, 2);              // diameter**2
    double S;                           // nematic order parameter
    double dl = 0.1;                    // maximum proposed displacement for each component
    double dn = 0.05;                   // maximum proposed change in direction for each component
    double drho = 0.001;                // maximum proposed change in rho
    double betaP = 1;                   // pressure P, temperature P/(k_B P T)
    double pvol = 0.5;                  // probability to propose a change in volume
    double rho;                         // density of the system (volume occupied by particles / total volume)
    double rho_proposed;                // density of the system (volume occupied by particles / total volume)
    double v0;                          // volume of one particle
    double V;                           // volume of the box
    double V_proposed;                  // proposed volume of the box
    double rho_cp;                      // close packed density 
    double accept_rate;                 // acceptance rate
    particle p_temp;                    // temporary particle holder
    vector<double> box_l;               // coordinates of box corners, one is at the origin
    vector<particle> part;              // main particle matrix
    string folder;                      // folder name to ouput coords

    simulation(string f, double l) : folder(f), L(l) 
    {
        v0 = M_PI * (L * pow(D, 2) / 4 + pow(D, 3) / 6);
        rho_cp = 2 / (sqrt(2) + (L/D)*sqrt(3));
    };

    inline void write_config(int step_no){write_config_to_file(step_no, N, box_l, part, L, D, folder);}

    // random integers from in the interval [0, N-1] 
    inline double rint(){return (int)(dsfmt_genrand() * N);}

    void fcc_config(int N_side, double rho_initial){
        rho = rho_initial;
        auto [part_c, box_l_c, N_c] = fcc_config_initial(N_side, L, rho_initial);
        part  = part_c;
        box_l = box_l_c;
        N     = N_c;
        V = box_l[0] * box_l[1] * box_l[2];
    }

    void update_S(){
        S = 0;
        vector<double> director(ndim);
        for(int i=0; i<N; i++){
            for(int j=0; j<ndim; j++){
                director[j] += part[i].dir[j];
            }
        }
        normalize(director);

        for(int i=0; i<N; i++){
            S += pow(v_sprod(director, part[i].dir), 2);
        }

        double Sb = (3*S / (double)N - 1)/2;
        cout << Sb << endl;
    }

    void scale(double scale_factor){
        for(int j=0; j<ndim; j++){
            for(int i=0; i<N; i++){
                part[i].pos[j] *= scale_factor;
            }
            box_l[j] *= scale_factor;
        }
    }

    double scale_factor(double r){
        V_proposed = N * v0 / (r * rho_cp);
        return cbrt(V_proposed / V);
    }

    void propose_NVT(){
        pmoved = rint();
        if(rdouble() < 0.5){
            for(int j=0; j<ndim; j++){
                p_temp.pos[j] = part[pmoved].pos[j] + (2 * rdouble() - 1) * dl;
                p_temp.pos[j] -= box_l[j] * floor(p_temp.pos[j]/box_l[j]);
            }
        }
        else{
            for(int j=0; j<ndim; j++){
                p_temp.dir[j] = part[pmoved].dir[j] + (2 * rdouble() - 1) * dn;
            }
            normalize(p_temp.dir);
        }
    }

    void NVT_step(){
        propose_NVT();
        metropolis_acceptance_NVT();
    }

    void NPT_step(){
        if(rdouble() < pvol){
            rho_proposed = rho + drho * (2 * rdouble() - 1);

            if(rho_proposed < 1e-4){
                return;
            }

            double factor = scale_factor(rho_proposed);
            scale(factor);

            /*
            cout << endl;
            cout << "betaP:     " << betaP << endl;
            cout << "V:         " << V << endl;
            cout << "delta_rho: " << rho_proposed - rho << endl;
            cout << "delta_V:   " << V_proposed - V << endl; 
            cout << "acc:       " << exp(-betaP*(V_proposed - V) + N * log(V_proposed/V)) << endl;
            cout << "acc1_term: " << -betaP*(V_proposed - V) << endl;
            cout << "acc2_term: " <<  N * log(V_proposed/V)  << endl;
            */

            if(!exists_general_overlap() && (rdouble() < exp(-betaP*(V_proposed - V) + N * log(V_proposed/V)))){
                //cout << "accepted!" << endl;
                V = V_proposed;
                rho = rho_proposed;
                accept++;
            }
            else{
                scale(1/factor);
            }
        }
        else{
            //cout << "NVT" << endl;
            NVT_step();
        }
    }

    bool they_overlap(particle &p1, particle &p2){
        return (dist_rods(p1, p2, box_l, L) < D2);
    }

    bool exists_general_overlap(){
        for(int i=0; i<N-1; i++){
            for(int j=i+1; j<N; j++){
                if(they_overlap(part[i], part[j]))
                    return true;
            }
        }
        return false;
    }

    bool exists_overlap(){
        for(int i=0; i<N; i++){
            if(i != pmoved){
                if(they_overlap(part[i], p_temp))
                    return true;
            }
        }
        return false;
    }

    void metropolis_acceptance_NVT(){
        if(!exists_overlap()){
            accept_proposal_NVT();
        }
    }

    void accept_proposal_NVT(){
        accept++;
        accept_rate = (double)accept/step;
        part[pmoved] = p_temp;
    }

};


