// prints 1D vector 
void print(std::vector<double> &v){
    for(int i=0; i<v.size(); i++){
        std::cout << v[i] << std::endl;
    }
}

// prints 2D vector 
void print(std::vector<std::vector<double>> &v){
    for(int i=0; i<v.size(); i++){
        std::cout << std::left;
        for(int j=0; j<v[0].size(); j++){
            std::cout << std::setw(20) << v[i][j];
        }
        std::cout << std::endl;
    }
}

// prints 3D vector 
void print(std::vector<std::vector<std::vector<double>>> &v){
    for(int i=0; i<v.size(); i++){
        std::cout << std::left;
        for(int j=0; j<v[0].size(); j++){
            for(int k=0; k<v[0][0].size(); k++){
                std::cout << std::setw(20) << v[i][j][k];
            }
        std::cout << std::endl;
        }
    std::cout << std::endl;
    }
}
