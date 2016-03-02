//
// Created by philipp on 11.01.16.
//

#ifndef HASHING_PERFECTHASHFUNCTION_H
#define HASHING_PERFECTHASHFUNCTION_H

#include "definitions.h"
#include "InputData.h"

class PerfectHashFunction {

private:
    // TODO smaller domains, stack instead of heap?, comments ignore arithmetic coding...
    // general parameters
    //ULLONG* _c;   // additional number of bits for 1-universal hash functions
    unsigned short _k;   // universe U=({0,1}^k)^l, needed for split of input x (especially if there are leading zero)
    unsigned short _l;   // needed for navigation in _h_split_coeffs
    // general information
    ULLONG _m;   // number of buckets, needed for navigation in mi, _h_mod_mask, ...
    ULLONG* _offset; // [m+1] -> [m_0+m_1+m_2+...+m_(m-1)+1], offset of bucket i, offset[m] first position after buckets
    // TODO offset[0] = 0 trivial, not necessary -> save as _offset_ and define function/... _offset? BUT computational overhead to test == 0
    //ULLONG* _mi;  // [0] .. [m-1], sizes of buckets, needed for computation of f^i_j, can be computed on demand: mi = offset[i+1]-offset[i]
    // h_split
    ULLONG _h_split_mod_mask;    // mask to compute MOD 2^(k+log(m)+c), value 2^(k+log(m)+c)-1
    // short _h_split_div_bits;    // k = number of bits for shift to compute DIV 2^k
    ULLONG* _h_split_coeffs;      // [l+1] -> [2^(k+log(m)+c)]=[_h_split_mod_mask + 1]
    // 1-universal hash functions h^i_j, 0<=i<m, 0<=j<2
    // ordering: h^0_0, h^0_1, h^1_0, h^1_1, h^2_0, h^2_1, ...
    ULLONG _h_mod_mask;      // mask to compute MOD 2^(k+log(r)+c), value 2^(k+log(r)+c)-1, same for all HFs
    //short _h_div_bits = _h_split_div_bits = k
    ULLONG _tab_rows;                // [_tab_rows] domain of HFs (=r)
    ULLONG* _h_coeffs;        // [(2*m)*(l+1)] -> [2^(k+log(_tab_rows)+c)]=[_h_mod_mask + 1]
    // tables T^0_0, T^0_1, T^1_0, T^1_1, T^2_0, T^2_1
    // ordering: tables = columns, row-wise representation
    //_tab_rows = number of rows of each table
    unsigned short _tab_width;       // = b = number of bits of each table entry (b = log(max ni) + c)
    ULLONG* _random_table;    // [6*_tab_rows] -> [2^(_tab_width)]
    ULLONG* _random_factor;   // si's, [m] (or [3*m]) -> [2^(_tab_width)]
    // selection arrays G_i where G_i: [m_i] -> [3], i in [m]
    // ordering: G_0[0], G_0[1], G_0[2], ..., G_0[m_0-1], G_1[0], G_1[1], ..., G_(m-1)[m_(m-1)-1]
    // G_i[j] = _g[offset[i]+j]
    unsigned char* _g;       // [offset[m]] -> [3]; 8 bit <-> 4 values
    bool _debug_mode;       // indicates whether to output debug notes or not

    void _configure(Configuration &config, ULLONG data_length); //step 1 + bit masks
    void _createUhf(ULLONG* coeffs, mt19937* rng, uniform_int_distribution<ULLONG>* dist); //step 2
    ULLONG _evalUhf(ULLONG key, ULLONG* coeff, ULLONG modMask, ULLONG modulus);
    bool _split(Configuration &config, InputData *data, InputData *bucket_data,
                ULLONG *&bucket_offsets, ULLONG *max_bucket_size, ULLONG *max_mi, Statistics &stats); //steps 4-7 + _tab_width
    void _createGoodPairs(Configuration &config, InputData *bucket_data, ULLONG *bucket_offsets, ULLONG max_bucket_size, mt19937* rng,
                          uniform_int_distribution<ULLONG>* dist, Statistics &stats); //step 8
    void _createRandomTables(mt19937* rng, uniform_int_distribution<ULLONG>* dist); //step 9
    void _createRandomFactor(ULLONG bucket_num, mt19937* rng, uniform_int_distribution<ULLONG>* dist); //step 10.1
    void _computeGij(ULLONG bucket_num, ULLONG *acyclicity_test_array, InputData *bucket_data, ULLONG bucket_offset, ULLONG bucket_size); //step 10.2, 11-12
    bool _isCyclic(ULLONG bucket_num, ULLONG *acyclicity_test_array, ULLONG bucket_size); //step 13+_g, reset acyclicity_test_array to zero...
    void _peelOf(ULLONG edge_index, ULLONG vertex_index, ULLONG *acyclicity_test_array,
                 ULLONG bucket_size, ULLONG *queue, ULLONG &next_queue_index, ULLONG *edgesOf,
                 ULLONG *cEdgesOf, ULLONG max_length, ULLONG mi, unsigned char *removed);
    void _clear(); // delete data
    void _computeSizes(Statistics &stats);

public:
    virtual ~PerfectHashFunction() { _clear(); }

public:
    PerfectHashFunction(Configuration &config, InputData *data, Statistics &stats);
    // TODO int instead of ULLONG? only one bit/one number needed to express "not found": ULLONG_MAX?
    ULLONG evaluate(ULLONG x);
    ULLONG getRange();
    ULLONG getSizeInBytes();
    ULLONG getCompactSizeInBytes();
};

#endif //HASHING_PERFECTHASHFUNCTION_H
