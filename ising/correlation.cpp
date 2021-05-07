#include <iostream>
#include "mt19937.h"

using namespace std;

int main(){
    dsfmt_seed(time(NULL));

    // size of the array
    int size = 1000;
    
    // creating a random array just to test the correlation function
    cout << "array:" << endl;
    double m[size];
    for(int i=0; i<size; i++){
        m[i] = 1;
        cout << m[i] << endl;
    }

    // this is the correlation calculation for a general_array "m" of length "size"
    int t_max = size - 1;
    double correlation[t_max]; 
    double sum_correlation = 0;
    double sum1, sum2, sum3, factor;

    cout << "\ncorrelation:" << endl;
    for(int t=0; t<t_max; t++){
        sum1 = 0;
        sum2 = 0;
        sum3 = 0;

        for(int s=0; s<t_max-t+1; s++){
            sum1 += m[s] * m[s+t];
            sum2 += m[s];
            sum3 += m[s+t];
        }

        factor = (double)(t_max - t);
        correlation[t] = (sum1 - sum2*sum3/factor) / factor;
        sum_correlation += correlation[t];
        cout << correlation[t] << endl;
    }

    cout << "\nestimate correlation time:" << endl;
    cout << sum_correlation / correlation[0] << endl;



}
