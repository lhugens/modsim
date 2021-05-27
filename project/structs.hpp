using namespace std;

#define ndim 3

struct particle{
    vector<double> pos;
    vector<double> dir;
    particle() : pos(ndim), dir(ndim) {}
};

struct simulation{
    double box_l;
    double L = 1;
    double D = 0.5;

    simulation(double BOX_L) : box_l(BOX_L)
    {}

    double dist_rods(particle & p1, particle & p2){
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
    
};
