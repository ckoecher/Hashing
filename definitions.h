/*
 * Hashing - A simple implementation of Split-And-Share-Hashing.
 * Copyright (C) 2016  Philipp Schlag, Chris Köcher
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
#include <sys/stat.h>

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
/**
 * This struct stores the configuration data (parameters etc.) for current perfect hash function creation process.
 */
struct Configuration {
    /**
     * The bit length of one key fragment.
     */
    unsigned short k;

    /**
     * The number of fragments that a key is split into.
     */
    unsigned short l;

    /**
     * The coefficient of m, which is the number of buckets and so the range of the split functions image.
     */
    double m_coeff;

    /**
     * The exponent of m, which is the number of buckets and so the range of the split functions image.
     * Let n be the number of keys. Then m = m_coeff * n ^ m_exp.
     */
    double m_exp;

    /**
     * The number of additional bits for computation of 1-universal hash functions to obtain (approx.) full randomness.
     */
    unsigned short additional_bits_uhf;

    /**
     * The maximal number of tries to split the keys into small buckets.
     */
    unsigned short num_of_tries_split;

    /**
     * The maximal number of tries to find a good pair of 1-universal hash functions (for each bucket).
     */
    unsigned short num_of_tries_goodpairs;

    /**
     * The coefficient of m_i, which is the range of the image of the buckets hash function.
     * Let n_i be the size of the bucket i. Then m_i = mi_coeff * n_i.
     */
    double mi_coeff;

    /**
     * The coefficient in the calculation of the number of rows in the tables T^i_j.
     */
    double tab_rows_coeff;

    /**
     * The exponent in the calculation of the number of rows in the tables T^i_j.
     * Let N be the maximal bucket size. Then the table will have tab_rows_coeff * (N ^ tab_rows_exp) rows.
     */
    double tab_rows_exp;

    /**
     * The number of additional bits for each table entry in T^i_j.
     */
    unsigned short additional_bits_tab;

    /**
     * The maximal number of tries to find tables T^i_j of random values to make the bucket hash functions perfect.
     */
    unsigned short num_of_tries_random_tab;

    /**
     * The maximal number of tries to find a random factor to make the bucket hash functions perfect (for each bucket).
     */
    unsigned short num_of_tries_random_si;

    /**
     * The seed for the random number generator. Consider that mt19937::result_type is the same as unsigned long.
     */
    mt19937::result_type seed;

    /**
     * Indicates wheth to output debugging output (true) or not (false).
     */
    bool debug_mode;
};

/**
 * This struct stores some statistics of current perfect hash function creation.
 */
struct Statistics {
    // Some general data
    /**
     * The number of clock ticks per second of the machine that is running this software.
     */
    ULLONG clocks_per_sec = CLOCKS_PER_SEC;

    /**
     * The number of keys in data file.
     */
    ULLONG num_of_keys = 0;

    /**
     * The range of the image of the perfect hash function.
     */
    ULLONG range_of_phf = 0;


    // The actual used bytes
    /**
     * The size of the whole PHF.
     */
    ULLONG size_in_bytes = 0;

    /**
     * The size of k (bit length per key fragment) and l (number of key fragments per key).
     */
    ULLONG size_in_bytes_general = 0;

    /**
     * The size of the coefficients of the split function.
     */
    ULLONG size_in_bytes_split_uhf = 0;

    /**
     * The size of the offset array.
     */
    ULLONG size_in_bytes_offsets = 0;

    /**
     * The size of the good hash function pairs h^i_j
     */
    ULLONG size_in_bytes_good_uhf_pairs = 0;

    /**
     * The size of PerfectHashFunction::_tab_width.
     */
    ULLONG size_in_bytes_random_width = 0;

    /**
     * The size of the random tables T^i_j.
     */
    ULLONG size_in_bytes_random_table = 0;

    /**
     * The size of the additional random factors s_i.
     */
    ULLONG size_in_bytes_random_factor = 0;

    /**
     * The size of the G_i arrays.
     */
    ULLONG size_in_bytes_g_array = 0;


    // The needed bytes
    /**
     * The size of the whole PHF.
     */
    ULLONG compact_size_in_bytes = 0;

    /**
     * The size of k (bit length per key fragment) and l (number of key fragments per key).
     */
    ULLONG compact_size_in_bytes_general = 0;

    /**
     * The size of the coefficients of the split function.
     */
    ULLONG compact_size_in_bytes_split_uhf = 0;

    /**
     * The size of the offset array.
     */
    ULLONG compact_size_in_bytes_offsets = 0;

    /**
     * The size of the good hash function pairs h^i_j
     */
    ULLONG compact_size_in_bytes_good_uhf_pairs = 0;

    /**
     * The size of PerfectHashFunction::_tab_width.
     */
    ULLONG compact_size_in_bytes_random_width = 0;

    /**
     * The size of the random tables T^i_j.
     */
    ULLONG compact_size_in_bytes_random_table = 0;

    /**
     * The size of the additional random factors s_i.
     */
    ULLONG compact_size_in_bytes_random_factor = 0;

    /**
     * The size of the G_i arrays.
     */
    ULLONG compact_size_in_bytes_g_array = 0;


    // The statistics of the whole creation process
    /**
     * The number of ticks until start of the creation process.
     */
    clock_t creation_start = 0;

    /**
     * The number of ticks until start of the creation process.
     */
    clock_t creation_end = 0;

    /**
     * The number of ticks needed for creation process.
     */
    clock_t creation_time = 0;

    /**
     * The number of ticks needed I/O stuff on whole creation process.
     */
    clock_t creation_io = 0;

    /**
     * Indicates whether the creation process was successful (true) or not (false).
     */
    bool creation_success = false;


    // The statistics of setup process
    /**
     * The number of ticks until start of the setup process.
     */
    clock_t setup_start = 0;

    /**
     * The number of ticks until end of the setup process.
     */
    clock_t setup_end = 0;

    /**
     * The number of ticks needed for setup process.
     */
    clock_t setup_time = 0;

    /**
     * The number of ticks needed I/O stuff on setup process.
     */
    clock_t setup_io = 0;

    /**
     * Indicates whether the setup process was successful (true) or not (false).
     */
    bool setup_success = false;


    // The statistics of splitting process
    /**
     * The number of ticks until start of the splitting process.
     */
    clock_t split_start = 0;

    /**
     * The number of ticks until end of the splitting process.
     */
    clock_t split_end = 0;

    /**
     * The number of ticks needed for splitting process.
     */
    clock_t split_time = 0;

    /**
     * The number of ticks needed I/O stuff on splitting process.
     */
    clock_t split_io = 0;

    /**
     * The number of tries needed to split the data.
     */
    unsigned short split_tries = 0;

    /**
     * Indicates whether the splitting process was successful (true) or not (false).
     */
    bool split_success = false;


    // The statistics on split function
    /**
     * The number of buckets.
     */
    ULLONG num_of_buckets = 0;

    /**
     * The maximal size of a bucket.
     */
    ULLONG max_bucket_size = 0;

    /**
     * The minimal size of a bucket.
     */
    ULLONG min_bucket_size = 0;

    /**
     * The average size of one bucket.
     */
    long double avg_bucket_size = 0.0;


    // The statistics of creation process of good function pairs h^i_j
    /**
     * The number of ticks until start of the creation process of the good function pairs.
     */
    clock_t goodpairs_start = 0;

    /**
     * The number of ticks until end of the creation process of the good function pairs.
     */
    clock_t goodpairs_end = 0;

    /**
     * The number of ticks needed for creation process of the good function pairs.
     */
    clock_t goodpairs_time = 0;

    /**
     * The number of ticks needed for I/O on creation process of the good function pairs.
     */
    clock_t goodpairs_io = 0;

    /**
     * The total number of tries needed to find good function pairs.
     */
    ULLONG goodpairs_total_tries = 0;

    /**
     * Indicates whether we've found good function pairs for each bucket (true) or not (false).
     */
    bool goodpairs_success = false;


    // The statistics for finding process of each buckets hash functions
    /**
     * The number of ticks until start of the finding process of each buckets hash functions.
     */
    clock_t buckets_start = 0;

    /**
     * The number of ticks until end of the finding process of each buckets hash functions.
     */
    clock_t buckets_end = 0;

    /**
     * The number of ticks needed for finding process of each buckets hash functions.
     */
    clock_t buckets_time = 0;

    /**
     * The number of ticks needed for I/O on finding process of each buckets hash function.
     */
    clock_t buckets_io = 0;

    /**
     * The number of tries to find random hash tables T^i_j.
     */
    unsigned short random_tab_tries = 0;

    /**
     * The total number of tries to find random factors s_i.
     */
    ULLONG random_si_total_tries = 0;

    /**
     * Indicates whether we've created the hash functions for each bucket (true) or not (false).
     */
    bool buckets_success = false;


    // The statistics of evaluation process
    /**
     * The number of ticks until start of the evaluation process.
     */
    clock_t eval_start = 0;

    /**
     * The number of ticks until end of the evaluation process.
     */
    clock_t eval_end = 0;

    /**
     * The number of ticks needed for evaluation process.
     */
    clock_t eval_time = 0;

    /**
     * The number of ticks needed for I/O on evaluation process.
     */
    clock_t eval_io = 0;

    /**
     * The average number of ticks needed for evaluation of hash function on one key.
     */
    long double avg_eval_time = 0.0;

    /**
     * The average number of ticks needed for I/O on evaluation of hash function on one key.
     */
    long double avg_eval_io = 0.0;

    /**
     * Indicates whether the perfect hash function is injective (true) or not (false).
     */
    bool eval_success = false;
};

#endif //HASHING_DEFINITIONS_H
