//
// Created by philipp on 11.01.16.
//

#ifndef HASHING_PERFECTHASHFUNCTION_H
#define HASHING_PERFECTHASHFUNCTION_H

#include "definitions.h"
#include "InputData.h"

class PerfectHashFunction {

private:
    // general information
    unsigned short _k;          // bit length of one key fragment
    unsigned short _l;          // needed for navigation in _h_split_coeffs
    // split information
    ULLONG _m;                  // number of buckets the set of keys will be split into
    ULLONG *_offset = nullptr;            // [m+1] -> [m_0+m_1+m_2+...+m_(_m-1)+1], offset of range of bucket i, offset[m] first position after whole range
    ULLONG _h_split_mod_mask;   // [_h_split_mod_mask + 1] intermediate range of 1-universal hash function h_split
    ULLONG *_h_split_coeffs = nullptr;    // [_l+1] -> [_h_split_mod_mask + 1], coefficients of 1-universal hash function h_split
    // good 1-universal hash function pairs for each bucket information
    // 1-universal hash functions h^i_j where 0<=i<_m, 0<=j<2
    // ordering: h^0_0, h^0_1, h^1_0, h^1_1, h^2_0, h^2_1, ...
    ULLONG _h_mod_mask;         // [_h_mod_mask + 1] intermediate range of good 1-universal hash function pairs, same for each bucket
    ULLONG _tab_rows;           // [_tab_rows] range of good 1-universal hash function pairs, same for each bucket
    ULLONG *_h_coeffs = nullptr;          // [2*_m*(_l+1)] -> [_h_mod_mask + 1], coefficients of good 1-universal hash function pairs
    // shared randomness information
    // tables T^0_0, T^0_1, T^1_0, T^1_1, T^2_0, T^2_1
    // ordering: tables = columns, line by line representation
    // _tab_rows = number of rows of each table
    unsigned short _tab_width;  // number of bits of each table entry
    ULLONG *_random_table = nullptr;      // [6*_tab_rows] -> [2^_tab_width], table with random values for perfect bucket hash functions
    ULLONG *_random_factor = nullptr;     // [_m] -> [2^_tab_width], array with random factors for perfect bucket hash functions (one per bucket)
    // selection arrays G_i where G_i: [m_i] -> [3], i in [_m]
    // ordering: G_0[0], G_0[1], G_0[2], ..., G_0[m_0-1], G_1[0], G_1[1], ..., G_(_m-1)[m_(_m-1)-1]
    // G_i[j] = _g[_offset[i]+j]
    unsigned char *_g = nullptr;          // [_offset[_m]] -> [3], selects one of three hash values to guarantee injectivity, 8 bits <-> 4 values
    // debug information
    bool _debug_mode;           // indicates whether to output debug notes or not

    void _configure(Configuration &config, ULLONG data_length); // step 1 + bit masks
    void _createUhf(ULLONG *coeffs, mt19937 *rng, uniform_int_distribution<ULLONG> *dist); // step 2
    ULLONG _evalUhf(ULLONG key, ULLONG *coeff, ULLONG modMask, ULLONG modulus);

    bool _split(Configuration &config, InputData *data, InputData *bucket_data,
                ULLONG *&bucket_offsets, ULLONG *max_bucket_size, ULLONG *max_mi,
                Statistics &stats); // steps 4-7 + _tab_width
    void _createGoodPairs(Configuration &config, InputData *bucket_data, ULLONG *bucket_offsets, ULLONG max_bucket_size,
                          mt19937 *rng,
                          uniform_int_distribution<ULLONG> *dist, Statistics &stats); // step 8
    void _createRandomTables(mt19937 *rng, uniform_int_distribution<ULLONG> *dist); // step 9
    void _createRandomFactor(ULLONG bucket_num, mt19937 *rng, uniform_int_distribution<ULLONG> *dist); // step 10.1
    void _computeGij(ULLONG bucket_num, ULLONG *acyclicity_test_array, InputData *bucket_data, ULLONG bucket_offset,
                     ULLONG bucket_size); // step 10.2, 11-12
    bool _isCyclic(ULLONG bucket_num, ULLONG *acyclicity_test_array,
                   ULLONG bucket_size); // step 13+_g, reset acyclicity_test_array to zero...
    void _peelOf(ULLONG edge_index, ULLONG vertex_index, ULLONG *acyclicity_test_array,
                 ULLONG bucket_size, ULLONG *queue, ULLONG &next_queue_index, ULLONG *edgesOf,
                 ULLONG *cEdgesOf, ULLONG max_length, ULLONG mi, unsigned char *removed);
    void _computeSizes(Statistics &stats);
    void _clear(); // delete data

public:
    virtual ~PerfectHashFunction() { _clear(); }

public:
    PerfectHashFunction(Configuration &config, InputData *data, Statistics &stats);

    ULLONG evaluate(ULLONG x);

    ULLONG getRange();

    ULLONG getSizeInBytes();

    ULLONG getCompactSizeInBytes();
};

#endif //HASHING_PERFECTHASHFUNCTION_H
