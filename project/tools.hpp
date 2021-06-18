struct particle{
    vector<double> pos;
    vector<double> dir;
    particle() : pos(ndim), dir(ndim) {}
};


inline double rdouble(){
    return dsfmt_genrand();
}


double v_sprod(vector<double> &v1, vector<double> &v2){
    double sum = 0;
    for(int i=0; i<v1.size(); i++){
        sum += v1[i]*v2[i];
    }
    return sum;
} 

double v_norm2(vector<double> &v){
    double sum = 0;
    for(int i=0; i<ndim; i++){
        sum += pow(v[i], 2);
    }
    return sum;
}

void normalize(vector<double> &v){
    double normm = sqrt(v_norm2(v));
    for(int i=0; i<ndim; i++){
        v[i] /= normm;
    }
}

void print(particle &p){
    cout << left;
    for(int i=0; i<ndim; i++){
        cout << setw(15) << p.pos[i];
    }
    for(int i=0; i<ndim; i++){
        cout << setw(15) << p.dir[i];
    }
    cout << endl;
}


void print(vector<particle> part){
    for(auto p : part){
        print(p);
    }
}

void print(vector<double> &v){
    for(auto x : v){
        cout << x << endl;
    }
}

void print(vector<vector<double>> &v){
    for(auto x : v){
        for(auto y : x){
            cout << left << setw(10) << y;
        }
        cout << endl;
    }
}


static inline double sign(double a, double b){ 
    a = fabs(a);
    return (b < 0) ? -a : a;
}

double dist_rods(particle &p1, particle &p2, vector<double> &box_l, double &L){
    vector<double> r12(ndim);

    for(int i=0; i<ndim; i++){
        r12[i] = p1.pos[i] - p2.pos[i];
    }

    vector<double> min_r12(ndim);
    double fx;

    for(int i=0; i<ndim; i++){
        fx = fabs(r12[i]);
        min_r12[i] = ( (fx < box_l[i] - fx)? r12[i] : (r12[i] - ((r12[i] > 0)? box_l[i] : -box_l[i])));
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
    double cos_the = n[2]; 
    double sin_the = sqrt(pow(n[0], 2) + pow(n[1], 2));
    bool gimbal_lock = sin_the < 1e-7;
    double cos_phi = gimbal_lock ? 1.0 : n[0] / sin_the;
    double sin_phi = gimbal_lock ? 0.0 : n[1] / sin_the;
    /*
    cout << "cos_the " << cos_the << endl; 
    cout << "sin_the " << sin_the << endl; 
    cout << "cos_phi " << cos_phi << endl; 
    cout << "sin_phi " << sin_phi << endl; 
    cout << endl;
    */

    vector<vector<double>> rot_matrix {{cos_phi*cos_the, -sin_phi, -n[0]},
                                       {sin_phi*cos_the,  cos_phi, -n[1]},
                                       {sin_the,                0,  n[2]}};
 
    //print(rot_matrix);
    return rot_matrix;
}

void write_config_to_file(int step, int &N, vector<double> &box_l, vector<particle> &part, double &L, double &D, string &folder){

    string step_str = to_string(step);
    string filename = folder + "_step_" + string(9 - step_str.length(), '0') + step_str + ".dat";

    ofstream file;
    file.open(filename);

    file << N << endl;

    for(int i=0; i<ndim; i++){
        file << left << setw(10) << 0.0 << setw(10) << box_l[i] << endl;
    }

    for(int i=0; i<part.size(); i++){
        vector<vector<double>> rot_matrix = rotation_matrix(part[i].dir);
        //cout << rot_matrix.size() << " " << rot_matrix[0].size() << endl;

        // particle type
        file << "SPHEROCYLINDER(" << L << ")  " << left;

        for(auto x : part[i].pos){
            file << setw(20) << x;
        }

        for(auto x : rot_matrix){
            for(auto y : x){
                file << setw(20) << y;
            }
        }

        file << endl;
    }
    file.close();
}

tuple<vector<particle>, vector<double>, int, double> file_config_initial(string filename){
    vector<particle> part;
    vector<double> box_l(ndim);
    particle p_temp;

    string trash;
    double g, n;
    int N;

    ifstream infile(filename);
    string line;

    infile >> N;

    // read box
    for (int i = 0; i < ndim; i++){
        infile >> g >> box_l[i];
    }

    for(int i=0; i<N; i++){
        infile >> trash;

        for(int j=0; j<ndim; j++){
            infile >> p_temp.pos[j];
        }
        infile >> g >> g >> n;
        p_temp.dir[0] = -n;
        infile >> g >> g >> n;
        p_temp.dir[1] = -n;
        infile >> g >> g >> n;
        p_temp.dir[2] = n;

    	part.push_back(p_temp);
    }

    unsigned first = trash.find("(");
    unsigned last = trash.find(")");
    string new_trash = trash.substr(first+1, last-first-1);
    double L = stod(new_trash);
    return {part, box_l, N, L};
}

tuple<vector<particle>, vector<double>, int> fcc_config_initial(int N_side, double L, double rho){

    vector<particle> part;
    vector<double> box_l(ndim);
    particle p_temp;
    double side = 2*L;

    p_temp.dir[2] = 1.0;
    box_l[0] = side*N_side + side;
    box_l[1] = side*N_side + side;
    box_l[2] = side*N_side + side;

    //generate cubic for 
    for(int i = 0; i < N_side + 1; i++){
        for(int j = 0; j < N_side + 1; j++){
            for(int l = 0; l < N_side + 1; l++){
                p_temp.pos[0] = (double)i*side;
                p_temp.pos[1] = (double)j*side;
                p_temp.pos[2] = (double)l*side;
    	        part.push_back(p_temp);
            }
        }   
    }
    //generate fcc
    for(int i = 0; i<N_side; i++){
        for(int j=0; j<N_side; j++){
            for(int l=0; l<N_side + 1; l++){
            	p_temp.pos[0] = (double)(i+0.5)*side;
    		p_temp.pos[1] = (double)(j+0.5)*side; 
    		p_temp.pos[2] = (double)l*side;
    		part.push_back(p_temp);
            }
        }
    }

    for(int i = 0; i<N_side; i++){
        for(int j=0; j<N_side + 1; j++){
            for(int l=0; l<N_side; l++){
            	p_temp.pos[0] = (double)(i+0.5)*side; 
    		p_temp.pos[1] = (double)j*side; 
    		p_temp.pos[2] = (double)(l+0.5)*side;
    		part.push_back(p_temp);
            }
        }
    }

    for(int i = 0; i<N_side + 1; i++){
        for(int j=0; j<N_side; j++){
            for(int l=0; l<N_side; l++){
            	p_temp.pos[0] = (double)i*side;
    		p_temp.pos[1] = (double)(j+ 0.5)*side; 
    		p_temp.pos[2] = (double)(l+0.5)*side;
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


