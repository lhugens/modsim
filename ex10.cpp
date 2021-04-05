#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdexcept>
#include <vector>

#include "mt19937.h"

using namespace std;

struct solver {
    vector<vector<int>> board;
    vector<vector<int>> var_entries;
    const int n_rows = 9;
    const int n_cols = 9;
    int n_var_entries;
    double beyta = 0;
    double energy;

    solver(vector<vector<int>> BOARD) {
        // init random seed
        size_t seed = time(NULL);
        dsfmt_seed(seed);

        // if input board is valid, initialize board and the matrix of 
        // variable entries (the ones with 0's, those will be our system essentially)  
        if ((BOARD.size() == n_rows) and (BOARD[0].size() == n_cols)) {
            board = BOARD;
            //print_board();
            for (int i = 0; i < n_rows; i++) {
                for (int j = 0; j < n_cols; j++) {
                    if (board[i][j] == 0) {
                        var_entries.push_back({ i, j });
                        board[i][j] = rand() * 9 + 1;
                    }
                }
            }
            n_var_entries = var_entries.size();
            get_energy();
        }
        else {
            throw invalid_argument("INVALID BOARD");
        }
    }

    // prints the current board state
    void print_board() {
        for (int i = 0; i < n_rows; i++) {
            cout << left;
            for (int j = 0; j < n_cols; j++) {
                cout << setw(3) << board[i][j];
            }
            cout << endl;
        }
    }

    inline double rand() {
        return dsfmt_genrand();
    }

    // prints any 1D int vector
    void print_vector(vector<int>& v) {
        cout << left;
        for (int i = 0; i < v.size(); i++)
            cout << setw(3) << v[i];
        cout << endl;
    }

    // prints any 2D int vector
    void print_matrix(vector<vector<int>>& matrix) {
        for (int i = 0; i < matrix.size(); i++) {
            cout << left;
            for (int j = 0; j < matrix[0].size(); j++) {
                cout << setw(3) << matrix[i][j];
            }
            cout << endl;
        }
    }

    // returns the number of duplicates in an int vector
    // assumes elemets satisfy 0 <= v[i] <= 9
    int duplicates(vector<int> v) {
        // the "+1" is to account for the "0" elements
        vector<int> histogram(n_rows + 1, 0);
        for (int i = 0; i < v.size(); i++)
            histogram[v[i]] += 1;
        int sum = 0;
        for (int i = 0; i < histogram.size(); i++) {
            if (histogram[i] > 1) {
                sum += histogram[i] - 1;
            }
        }
        return sum;
    }

    int duplicates_in_row(int i_row) {
        return duplicates(board[i_row]);
    }

    int duplicates_in_col(int j_col) {
        vector<int> temp(9, 0);
        for (int i = 0; i < n_rows; i++)
            temp[i] = board[i][j_col];
        return duplicates(temp);
    }

    int duplicates_in_block(int i_corner, int j_corner) {
        vector<int> temp(9, 0);
        int iter = 0;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                temp[iter] = board[i_corner + i][j_corner + j];
                iter++;
            }
        }
        return duplicates(temp);
    }

    // returns the "energy" of the board, which we define
    // to be the sum of duplicate elements in all the blocks,
    // rows and cols
    void get_energy() {
        energy = 0;

        // duplicates in all rows
        for (int i = 0; i < n_rows; i++)
            energy += duplicates_in_row(i);

        // duplicates in all columns
        for (int j = 0; j < n_cols; j++)
            energy += duplicates_in_col(j);

        // duplicates in all blocks
        // iterate through corners of the blocks
        for (int i_corner = 0; i_corner < 9; i_corner += 3) {
            for (int j_corner = 0; j_corner < 9; j_corner += 3)
                energy += duplicates_in_block(i_corner, j_corner);
        }
    }

    int change_rand_entry() {
        int entry_index = n_var_entries * rand();
        vector<int> entry = var_entries[entry_index];
        int i = entry[0];
        int j = entry[1];
        int i_corner = i - (i % 3);
        int j_corner = j - (j % 3);

        // determine energy contribution of an element before change 
        int old_entry_energy = 0;
        old_entry_energy += duplicates_in_row(i);
        old_entry_energy += duplicates_in_col(j);
        old_entry_energy += duplicates_in_block(i_corner, j_corner);

        // save old entry and generate a new random one
        int old_value = board[i][j];
        board[i][j] = rand() * 9 + 1;

        // determine energy contribution of an element after change 
        int new_entry_energy = 0;
        new_entry_energy += duplicates_in_row(i);
        new_entry_energy += duplicates_in_col(j);
        new_entry_energy += duplicates_in_block(i_corner, j_corner);

        double dE = new_entry_energy - old_entry_energy;

        if (rand() < exp(-beyta * dE)) {
            // accept
            energy += dE;
            return 1;
        }
        else {
            // reject
            board[i][j] = old_value;
            beyta = 1; //allow "uphill" to prevent beiing stuck
            return 0;
        }
    }
};

int main() {
    vector<vector<int>> board{
                            {7, 8, 0, 4, 0, 0, 1, 2, 0},
                            {6, 0, 0, 0, 7, 5, 0, 0, 9},
                            {0, 0, 0, 6, 0, 1, 0, 7, 8},
                            {0, 0, 7, 0, 4, 0, 2, 6, 0},
                            {0, 0, 1, 0, 5, 0, 9, 3, 0},
                            {9, 0, 4, 0, 6, 0, 0, 0, 5},
                            {0, 7, 0, 3, 0, 0, 0, 1, 2},
                            {1, 2, 0, 0, 0, 7, 4, 0, 0},
                            {0, 4, 9, 2, 0, 6, 0, 0, 7}
    };

    try {
        solver sol(board);
       
        cout << "initial energy " << sol.energy << endl;
       /* int accepted = 0;
        for(int i=0; i<10; i++)
            accepted += sol.change_rand_entry();
        cout << accepted << endl;
        cout << "final energy " << sol.energy << endl;*/
        double n = 1;
        char buffer[100];
        sprintf(buffer, "energy_1.dat");
        FILE* fp = fopen(buffer, "w");
        int accepted = 0;
        while (sol.energy != 0) {
            cout << sol.energy << endl;
            double e = sol.energy;
            fprintf(fp, "%lf\n", e);
            for (int i = 0; i < 100; i++) {
                accepted += sol.change_rand_entry();
                sol.beyta += i;
            }
        }
        sol.print_board();
        fclose(fp);
    }
    catch (invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}