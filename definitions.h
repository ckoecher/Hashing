#ifndef HASHING_DEFINITIONS_H
#define HASHING_DEFINITIONS_H

#include <iostream>     // cout, cin
#include <fstream>      // ifstream, ofstream
//#include <cstdint>      // uint_fast64_t
#include <climits>      // ULLONG_MAX
#include <string>       // string, stoull (to unsigned long long), stod (to double)
#include <random>       // Mersenne twister mt19937
#include "./inih/cpp/INIReader.h"
//#include <functional>   // bind
#include <assert.h>
#include <ctime> // clock_t
using namespace std;

//typedef uint_fast64_t UINT;
typedef unsigned long long int ULLONG;

#define ARR(array, rows, cols, row, col) array[row*cols+col]

// array of type unsigned char*
#define GETBITPAIR(array, tab, index) ((array[index >> 1] >> ( 6 - ((index & 1)<<2) - ((tab & 1)<<1) ) ) % 4)
inline unsigned char getBitPair(unsigned char* array, unsigned short tab, ULLONG index) {
    return (array[index >> 1] >> ( 6 - ((index & 1)<<2) - ((tab & 1)<<1) ) ) & 3;
}
#define INCBITPAIR(array, tab, index) (array[index >> 1] = array[index >> 1] + (1 << ((1-tab)<<1) << ((1-(index&1))<<2)))
#define ZEROBITPAIRS(array, tab, index) (array[index >> 1] = 0)
// TODO ATTENTION: does not only reset array_tab[index] to zero BUT complete char array element!

// NOT USABLE
#define CHARBITPAIR(array, index) ((array[index >> 2] >> 2*(3 - (index & 3))) % 4)
#define CHARBITPAIRTABS(array, tab, index) ((unsigned char)(array[index >> 1] >> ( 6 - ((index & 1)<<2) - ((tab & 1)<<1) ) ) % 4)
// last % 4 not replaceable by & 3 ???
// CHARBITPAIRTABS fastest solution

// slower...
#define CHARBITPAIR3(array, index) ((array[index/4] >> 2*(3 - (index % 4))) % 4)
#define CHARBITPAIR2(array, index) ((array[index >> 2] >> 2*(3 - (index % 4))) % 4)

// gets and sets a single bit in an array of chars (used as bitmaps)
#define GETBIT(array, index) (((array[index >> 3]) >> (index & 7)) & 1)
#define SETBIT(array, index, value) (array[index >> 3] ^= (-value ^ array[index >> 3]) & (1 << (index & 7)))

// gets and sets two bits in an array of chars (used as bitmap)
#define GETCHARBITPAIR(array, index) ((array[index >> 2]) >> ((index & 3) << 1) & 3)
//#define SETCHARBITPAIR(array, index, value) (array[index >> 2] ^= (-value ^ array[index >> 2]) & (3 << ((index & 3) << 1)))
#define SETCHARBITPAIR(array, index, value) (array[index>>2] = array[index>>2] & (255-(3<<((index&3)<<1))) | (value<<((index&3)<<1)))

struct Configuration {
    unsigned short k; // U = [(2^k)^l], k*l(<)=64, k>=1, l>=1
    unsigned short l;
    double m_coeff; // m = ceil(m_coeff * n^(m_exp))
    double m_exp;
    unsigned short additional_bits_uhf; // for 1-universal hash function
    unsigned short num_of_tries_split;
    unsigned short num_of_tries_goodpairs;
    double mi_coeff; // mi = ceil(c * ni)
    double tab_rows_coeff; // rows of tables T_i^j = tab_rows_coeff * n^(tab_rows_exp)
    double tab_rows_exp;
    unsigned short additional_bits_tab; // for table entries
    unsigned short num_of_tries_random_tab;
    unsigned short num_of_tries_random_si;
    //bool multiple_seeds; // use just one or more seeds for randomization of Mersenne Twister
    mt19937::result_type seed; // mt19937::result_type = unsigned long (= uint_fast32_t?)
    bool debug_mode; //enable or disable debugger output
};

struct Statistics {
    // TODO ..._time without IO?
    // TODO split = split_io + split_work?

    ULLONG clocks_per_sec = CLOCKS_PER_SEC;
    ULLONG num_of_keys = 0;
    ULLONG range_of_phf = 0;
    ULLONG size_in_bytes = 0;

    clock_t creation_start = 0;
    clock_t creation_end = 0;
    clock_t creation_time = 0;
    bool creation_success = false;

    clock_t setup_start = 0;
    clock_t setup_end = 0;
    clock_t setup_time = 0;
    bool setup_succuess = false;

    clock_t split_start = 0;
    clock_t split_end = 0;
    clock_t split_time = 0;
    unsigned short split_tries = 0;
    bool split_success = false;

    ULLONG num_of_buckets = 0;
    ULLONG max_bucket_size = 0;
    ULLONG min_bucket_size = 0;
    long double avg_bucket_size = 0.0;

    clock_t goodpairs_start = 0;
    clock_t goodpairs_end = 0;
    clock_t goodpairs_time = 0;
    ULLONG goodpairs_total_tries = 0;
    bool goodpairs_success = false;

    clock_t buckets_start = 0;
    clock_t buckets_end = 0;
    clock_t buckets_time = 0;
    unsigned short random_tab_tries = 0;
    ULLONG random_si_total_tries = 0;
    bool buckets_success = false;

    clock_t eval_start = 0;
    clock_t eval_end = 0;
    clock_t eval_time = 0;
    long double avg_eval_time = 0.0;
    bool eval_success = false;
};

//std::mt19937 generator (123);
//std::uniform_real_distribution<double> dis(0.0, 1.0);
//
//double randomRealBetweenZeroAndOne = dis(generator);

#endif //HASHING_DEFINITIONS_H
