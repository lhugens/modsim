struct particle{
    vector<double> pos;
    vector<double> dir;
    particle() : pos(ndim), dir(ndim) {}
};


double v_sprod(vector<double> &v1, vector<double> &v2){
    double sum = 0;
    for(int i=0; i<v1.size(); i++){
        sum += v1[i]*v2[i];
    }
    return sum;
} 

double v_norm2(vector<double> &v){
    double sum = 0;
    for(auto x : v){
        sum += pow(x, 2);
    }
    return sum;
}

void print(vector<double> &v){
    for(auto x : v){
        cout << x << endl;
    }
}

static inline double sign(double a, double b){ 
    a = fabs(a);
    return (b < 0) ? -a : a;
}

double dist_rods(particle &p1, particle &p2, double &box_l, double &L){
    vector<double> r12(ndim);

    for(int i=0; i<ndim; i++){
        r12[i] = p1.pos[i] - p2.pos[i];
    }

    vector<double> min_r12(ndim);
    double fx;

    for(int i=0; i<ndim; i++){
        fx = fabs(r12[i]);
        min_r12[i] = ( (fx < box_l - fx)? r12[i] : (r12[i] - ((r12[i] > 0)? box_l : -box_l)));
    }

    double rr   = v_norm2(min_r12);
    double rw1  = v_sprod(min_r12, p1.dir);
    double rw2  = v_sprod(min_r12, p2.dir);
    double w1w2 = v_sprod(p1.dir, p2.dir);
    double cc   = 1 - pow(w1w2, 2); 
    double xla, xmu;

    if(cc < 1e-6){
        if(rw1 && rw2){
            xla =  rw1 / 2;
            xmu = -rw2 / 2;
        }
        else
            return rr;
    } 
    else{
        xla = ( rw1 - w1w2*rw2) / cc;
        xmu = (-rw2 + w1w2*rw1) / cc;
    }

    double lh = L / 2;

    if (fabs(xla) > lh || fabs(xmu) > lh){
        if (fabs(xla) - lh > fabs(xmu) - lh){
            xla = sign(lh, xla);
            xmu = xla*w1w2 - rw2;
            if(fabs(xmu) > lh) 
                xmu = sign(lh, xmu);
            }
        else{
            xmu = sign(lh,xmu);
            xla = xmu*w1w2 + rw1;
            if (fabs(xla) > lh) 
                xla = sign(lh, xla);
        }
    }

    return rr + pow(xla, 2) + pow(xmu, 2) + 2*(xmu*rw2 - xla*(rw1 + xmu*w1w2));
}

vector<vector<double>> rotation_matrix(vector<double> &n){
    double cos_phi = sqrt(pow(n[0], 2) + pow(n[1], 2));
    double sin_phi = n[2]; 
    double cos_the = cos_phi ? n[0] / cos_phi : 1.0;
    double sin_the = cos_phi ? n[1] / cos_phi : 0.0;

    vector<vector<double>> rot_matrix {{n[0], -sin_the, -sin_phi*cos_the},
                                                 {n[1],  cos_the, -sin_phi*sin_the},
                                                 {n[2],        0,           cos_phi}};

    return rot_matrix;
}

void init_particle_file(ofstream &file, int N, double box_l){
    file << N << endl;

    for(int i=0; i<ndim; i++){
        file << left << setw(10) << 0.0 << setw(10) << box_l << endl;
    }
}

void write_particle_to_file(ofstream &file, particle &p, double &L, double &D){
    vector<vector<double>> rot_matrix = rotation_matrix(p.dir);

    // particle type
    file << "SPHEROCYLINDER(" << L << ")  " << left;

    for(auto x : p.pos){
        file << setw(10) << x;
    }

    for(auto x : rot_matrix){
        for(auto y : x){
            file << setw(10) << D*y;
        }
    }

    file << endl;
}







