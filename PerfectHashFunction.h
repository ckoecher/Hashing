//
// Created by philipp on 11.01.16.
//

#ifndef HASHING_PERFECTHASHFUNCTION_H
#define HASHING_PERFECTHASHFUNCTION_H

#include "definitions.h"

class PerfectHashFunction {

private:
    // TODO smaller domains, stack instead of heap?, comments ignore arithmetic coding...
    // general parameters
    //ULLONG* _c;   // additional number of bits for 1-universal hash functions
    short _k;   // universe U=({0,1}^k)^l, needed for split of input x (especially if there are leading zero)
    short _l;   // needed for navigation in _h_split_coeffs
    // general information
    ULLONG _m;   // number of buckets, needed for navigation in mi, _h_mod_mask, ...
    ULLONG* _offset; // [m] -> [m_0+m_1+m_2+...+m_(m-1)+1], offset of bucket i, offset[m] first position after buckets
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
    ULLONG _tab_rows;                // [_tab_rows] domain of HFs
    ULLONG* _h_coeffs;        // [(2*m)*(l+1)] -> [2^(k+log(_tab_rows)+c)]=[_h_mod_mask + 1]
    // tables T^0_0, T^1_0, T^2_0, T^0_1, T^1_1, T^2_2
    // ordering: tables = columns, row-wise representation
    //_tab_rows = number of rows of each table
    short _tab_width;       // = b = number of bits of each table entry (b = log(max ni) + c)
    ULLONG* _random_table;    // [6*_tab_rows] -> [2^(_tab_width)]
    ULLONG* _random_factor;   // si's, [m] (or [3*m]) -> [2^(_tab_width)]
    // selection arrays G_i where G_i: [m_i] -> [3], i in [m]
    // ordering: G_0[0], G_0[1], G_0[2], ..., G_0[m_0-1], G_1[0], G_1[1], ..., G_(m-1)[m_(m-1)-1]
    // G_i[j] = _g[offset[i]+j]
    unsigned char* _g;       // [offset[m]] -> [3]; 8 bit <-> 4 values

    void _configure(Configuration, ULLONG); //step 1 + bit masks
    void _createUhf(ULLONG*, mt19937*, uniform_int_distribution<ULLONG>*); //step 2
    ULLONG _evalUhf(ULLONG*, ULLONG);
    bool _split(Configuration, ULLONG*, ULLONG**, ULLONG*, ULLONG*, ULLONG*); //steps 4-7 + _tab_width
    void _createGoodPairs(ULLONG**, ULLONG*, mt19937*, uniform_int_distribution<ULLONG>*); //step 8
    void _createRandomTables(ULLONG, mt19937*, uniform_int_distribution<ULLONG>*); //step 9
    void _createRandomFactor(ULLONG, mt19937*, uniform_int_distribution<ULLONG>*); //step 10.1
    void _computeFij(ULLONG, ULLONG*, ULLONG, ULLONG*); //step 10.2
    void _computeGij(ULLONG*); //step 11-12
    bool _isCyclic(ULLONG, ULLONG*); //step 13+_g

public:
    virtual ~PerfectHashFunction() { }

public:
    PerfectHashFunction(Configuration config, ULLONG data_length, ULLONG *data);
    // TODO int instead of ULLONG? only one bit/one number needed to express "not found": ULLONG_MAX?
    ULLONG evaluate(ULLONG x);
};

#endif //HASHING_PERFECTHASHFUNCTION_H
