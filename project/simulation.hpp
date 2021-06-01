struct simulation{
    int N;                              // number of particles
    double L  = 1;                      // length of the spherocylinders
    double D  = 0.5;                    // Diameter of the cylinders
    double D2 = pow(D, 2);              // Diameter**2
    double dl = 0.1;                    // maximum proposed displacement for each component
    double dn = 0.1;                    // maximum proposed change in direction for each component
    vector<double> box_l;               // coordinates of box corners, one is at the origin
    vector<particle> part;              // main particle matrix
    vector<particle> part_proposed;     // main particle matrix

    simulation(int n) : N(n), part(n), part_proposed(n), box_l(ndim)
    {}

    inline void write_config(int step){write_config_to_file(step, N, box_l, part, L, D);}

    void propose(){
        for(int i=0; i<N; i++){
            for(int j=0; j<ndim; j++){
                part_proposed[i].pos[j] = part[i].pos[j] + rdouble() * dl;
                part_proposed[i].dir[j] = part[i].dir[j] + rdouble() * dn;
            }
            normalize(part_proposed[i].dir);
        }
    }

    void accept_proposal(){
        for(int i=0; i<N; i++){
            part[i] = part_proposed[i];
        }
    }

    bool they_overlap(particle &p1, particle &p2){
        return (dist_rods(p1, p2, box_l, L) < D2);
    }

    bool exists_overlap(){
        for(int i=0; i<N-1; i++){
            for(int j=i+1; j<N; j++){
                if(they_overlap(part_proposed[i], part_proposed[j]))
                    return true;
            }
        }
        return false;
    }

    void metropolis_acceptance(){
        if(!exists_overlap())
            accept_proposal();
    }
};


