#ifndef HASHING_DEFINITIONS_H
#define HASHING_DEFINITIONS_H

#include <iostream>     // cout, cin
#include <fstream>      // ifstream, ofstream
//#include <cstdint>      // uint_fast64_t
#include <climits>      // ULLONG_MAX
#include <string>       // string, stoull (to unsigned long long), stod (to double)
#include <random>       // Mersenne twister mt19937
//#include <functional>   // bind
#include <assert.h>
using namespace std;

//typedef uint_fast64_t UINT;
typedef unsigned long long int ULLONG;

#define ARR(array, rows, cols, row, col) array[row*cols+col]

struct Configuration {
    short k; // U = [(2^k)^l], k*l(<)=64, k>=1, l>=1
    short l;
    double m_coeff; // m = ceil(m_coeff * n^(m_exp))
    double m_exp;
    short additional_bits_uhf; // for 1-universal hash function
    double mi_coeff; // mi = ceil(c * ni)
    double tab_rows_coeff; // rows of tables T_i^j = tab_rows_coeff * n^(tab_rows_exp)
    double tab_rows_exp;
    short additional_bits_tab; // for table entries
    short num_of_tries_random_tab;
    short num_of_tries_random_si;
    //bool multiple_seeds; // use just one or more seeds for randomization of Mersenne Twister
    mt19937::result_type seed; // mt19937::result_type = unsigned long (= uint_fast32_t?)
};

//std::mt19937 generator (123);
//std::uniform_real_distribution<double> dis(0.0, 1.0);
//
//double randomRealBetweenZeroAndOne = dis(generator);

#endif //HASHING_DEFINITIONS_H
