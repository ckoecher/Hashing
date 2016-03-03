#ifndef HASHING_DEFINITIONS_H
#define HASHING_DEFINITIONS_H

#include <iostream>     // cout, cin
#include <fstream>      // ifstream, ofstream, fstream
#include <climits>      // ULLONG_MAX
#include <string>       // string, stoull (to unsigned long long), stod (to double)
#include <random>       // Mersenne twister mt19937
#include "./inih/cpp/INIReader.h"
#include <assert.h>
#include <ctime>        // clock_t

using namespace std;

// data type of keys
typedef unsigned long long int ULLONG;

// allows two-dimensional line by line navigation in an one-dimensional array
#define ARR(array, rows, cols, row, col) array[row*cols+col]

// gets and sets two bits in an array of unsigned chars (used as bitmap)
#define GETTABBITPAIR(array, tab, index) ((array[index >> 1] >> ( 6 - ((index & 1)<<2) - ((tab & 1)<<1) ) ) % 4)
#define INCTABBITPAIR(array, tab, index) (array[index >> 1] = array[index >> 1] + (1 << ((1-tab)<<1) << ((1-(index&1))<<2)))
#define ZEROTABBITPAIRS(array, index) (array[index >> 1] = 0)

// gets and sets a single bit in an array of unsigned chars (used as bitmap)
#define GETBIT(array, index) (((array[index >> 3]) >> (index & 7)) & 1)
#define SETBIT(array, index, value) (array[index >> 3] ^= (-value ^ array[index >> 3]) & (1 << (index & 7)))

// gets and sets two bits in an array of unsigned chars (used as bitmap)
#define GETCHARBITPAIR(array, index) ((array[index >> 2]) >> ((index & 3) << 1) & 3)
#define SETCHARBITPAIR(array, index, value) (array[index>>2] = array[index>>2] & (255-(3<<((index&3)<<1))) | (value<<((index&3)<<1)))

// stores configuration (parameters etc.) for current perfect hash function creation
struct Configuration {
    unsigned short k;                       // bit length of one key fragment
    unsigned short l;                       // number of fragments that a key is split into
    double m_coeff;                         // let n be the number of keys; then the set of keys will be
    double m_exp;                           // split into m_coeff * (n ^ m_exp buckets)
    unsigned short additional_bits_uhf;     // additional bits for computation of 1-universal hash functions to obtain (approx.) full randomness
    unsigned short num_of_tries_split;      // maximal number of tries to split the keys into small buckets
    unsigned short num_of_tries_goodpairs;  // maximal number of tries to find a good pair of 1-universal hash functions (for each bucket)
    double mi_coeff;                        // the keys in a bucket of size ni shall be hashed into a range of mi_coeff * ni values
    double tab_rows_coeff;                  // let N be the maximal bucket size; then the table of random values     rows of tables T_i^j = tab_rows_coeff * n^(tab_rows_exp)
    double tab_rows_exp;                    // will have tab_rows_coeff * (N ^ tab_rows_exp) rows
    unsigned short additional_bits_tab;     // additional bits for each table entry
    unsigned short num_of_tries_random_tab; // maximal number of tries to find a table of random values to make the bucket hash functions perfect
    unsigned short num_of_tries_random_si;  // maximal number of tries to find a random factor to make the bucket hash functions perfect (for each bucket)
    mt19937::result_type seed;              // seed for the random number generator; mt19937::result_type == unsigned long
    bool debug_mode;                        // enable or disable (debugger) output
};

// stores statistics of current perfect hash function creation
struct Statistics {
    ULLONG clocks_per_sec = CLOCKS_PER_SEC;
    ULLONG num_of_keys = 0;
    ULLONG range_of_phf = 0;

    ULLONG size_in_bytes = 0;
    ULLONG size_in_bytes_general = 0;
    ULLONG size_in_bytes_split_uhf = 0;
    ULLONG size_in_bytes_offsets = 0;
    ULLONG size_in_bytes_good_uhf_pairs = 0;
    ULLONG size_in_bytes_random_width = 0;
    ULLONG size_in_bytes_random_table = 0;
    ULLONG size_in_bytes_random_factor = 0;
    ULLONG size_in_bytes_g_array = 0;

    ULLONG compact_size_in_bytes = 0;
    ULLONG compact_size_in_bytes_general = 0;
    ULLONG compact_size_in_bytes_split_uhf = 0;
    ULLONG compact_size_in_bytes_offsets = 0;
    ULLONG compact_size_in_bytes_good_uhf_pairs = 0;
    ULLONG compact_size_in_bytes_random_width = 0;
    ULLONG compact_size_in_bytes_random_table = 0;
    ULLONG compact_size_in_bytes_random_factor = 0;
    ULLONG compact_size_in_bytes_g_array = 0;

    clock_t creation_start = 0;
    clock_t creation_end = 0;
    clock_t creation_time = 0;
    clock_t creation_io = 0;
    bool creation_success = false;

    clock_t setup_start = 0;
    clock_t setup_end = 0;
    clock_t setup_time = 0;
    clock_t setup_io = 0;
    bool setup_succuess = false;

    clock_t split_start = 0;
    clock_t split_end = 0;
    clock_t split_time = 0;
    clock_t split_io = 0;
    unsigned short split_tries = 0;
    bool split_success = false;

    ULLONG num_of_buckets = 0;
    ULLONG max_bucket_size = 0;
    ULLONG min_bucket_size = 0;
    long double avg_bucket_size = 0.0;

    clock_t goodpairs_start = 0;
    clock_t goodpairs_end = 0;
    clock_t goodpairs_time = 0;
    clock_t goodpairs_io = 0;
    ULLONG goodpairs_total_tries = 0;
    bool goodpairs_success = false;

    clock_t buckets_start = 0;
    clock_t buckets_end = 0;
    clock_t buckets_time = 0;
    clock_t buckets_io = 0;
    unsigned short random_tab_tries = 0;
    ULLONG random_si_total_tries = 0;
    bool buckets_success = false;

    clock_t eval_start = 0;
    clock_t eval_end = 0;
    clock_t eval_time = 0;
    clock_t eval_io = 0;
    long double avg_eval_time = 0.0;
    long double avg_eval_io = 0.0;
    bool eval_success = false;
};

#endif //HASHING_DEFINITIONS_H
