void print(std::vector<std::vector<double>> &v){
    for(int i=0; i<v.size(); i++){
        std::cout << std::left;
        for(int j=0; j<v[0].size(); j++){
            std::cout << std::setw(20) << v[i][j];
        }
        std::cout << std::endl;
    }
}

void print(std::vector<double> &v){
    for(int i=0; i<v.size(); i++){
        std::cout << v[i] << std::endl;
    }
}
