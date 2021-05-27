double v_sprod(std::vector<double> &v1, std::vector<double> &v2){
    double sum = 0;
    for(int i=0; i<v1.size(); i++){
        sum += v1[i]*v2[i];
    }
    return sum;
} 

double v_norm2(std::vector<double> &v){
    double sum = 0;
    for(auto x : v){
        sum += pow(x, 2);
    }
    return sum;
}

std::vector<double> v_subtract(std::vector<double> &v1, std::vector<double> &v2){
    std::vector<double> result(v1.size());
    for(int i=0; i<result.size(); i++){
        result[i] = v1[i] - v2[i];
    }
    return result;
}


void print(std::vector<double> &v){
    for(auto x : v){
        std::cout << x << std::endl;
    }
}

static inline double sign(double a, double b){ 
    a = fabs(a);
    return (b < 0) ? -a : a;
}
