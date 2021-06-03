struct simulation{
    int N;                              // number of particles
    double L  = 1;                      // length of the spherocylinders
    double D  = 0.5;                    // Diameter of the cylinders
    double D2 = pow(D, 2);              // Diameter**2
    double dl = 0.1;                    // maximum proposed displacement for each component
    double dn = 0.1;                    // maximum proposed change in direction for each component
    double density;                     // density of the system (volume occupied by particles / total volume)
    vector<double> box_l;               // coordinates of box corners, one is at the origin
    vector<particle> part;              // main particle matrix
    vector<particle> part_proposed;     // main particle matrix

    simulation(int n) : N(n), part(n), part_proposed(n), box_l(ndim)
    {}

    inline void write_config(int step){write_config_to_file(step, N, box_l, part, L, D);}

    void set_density(double dens){
        density = dens;

        double v0 = M_PI * (L * pow(D, 2) / 4 + pow(D, 3) / 6);
        double current_volume = box_l[0] * box_l[1] * box_l[2];
        double desired_volume = N * v0 / (2 * density/(sqrt(2) + (L/D)*sqrt(3)));
        double scale_factor   = cbrt(desired_volume / current_volume);

        for(int j=0; j<ndim; j++){
            for(int i=0; i<N; i++){
                part[i].pos[j] *= scale_factor;
            }
            box_l[j] *= scale_factor;
        }
    }

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

    void fcc_config(){

	int N_side = 4;
	int N_total = 0;

        vector<particle> part_fcc;

        box_l[0] = 2*L*N_side;
        box_l[1] = 2*L*N_side;
        box_l[2] = 2*L*N_side;

        particle p_temp;
        p_temp.dir[2] = 1.0;

        //generate cubic for 
        for(int i = 0; i < N_side + 1; i++){
            for(int j = 0; j < N_side + 1; j++){
                for(int l = 0; l < N_side + 1; l++){
		    N_total++;
                    p_temp.pos[0] = (double)i*2*L;
                    p_temp.pos[1] = (double)j*2*L;
                    p_temp.pos[2] = (double)l*2*L;
        	    part_fcc.push_back(p_temp);
                }
            }   
        }
        //generate fcc
        for(int i = 0; i<N_side; i++){
            for(int j=0; j<N_side; j++){
                for(int l=0; l<N_side + 1; l++){
			N_total += 1;
                	p_temp.pos[0] = (double)(i+0.5)*2*L;
			p_temp.pos[1] = (double)(j+0.5)*2*L; 
			p_temp.pos[2] = (double)l*2*L;
        		part_fcc.push_back(p_temp);
                }
            }
        }

        for(int i = 0; i<N_side; i++){
            for(int j=0; j<N_side + 1; j++){
                for(int l=0; l<N_side; l++){
			N_total += 1;
                	p_temp.pos[0] = (double)(i+0.5)*2*L; 
			p_temp.pos[1] = (double)j*2*L; 
			p_temp.pos[2] = (double)(l+0.5)*2*L;
        		part_fcc.push_back(p_temp);
                }
            }
        }

        for(int i = 0; i<N_side + 1; i++){
            for(int j=0; j<N_side; j++){
                for(int l=0; l<N_side; l++){
			N_total += 1;
                	p_temp.pos[0] = (double)i*2*L;
			p_temp.pos[1] = (double)(j+ 0.5)*2*L; 
			p_temp.pos[2] = (double)(l+0.5)*2*L;
        		part_fcc.push_back(p_temp);
                }
            }
        }

	write_config_to_file(0, N_total, box_l, part_fcc, L, D);

    }
};


