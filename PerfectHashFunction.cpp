//
// Created by philipp on 11.01.16.
//

#include "definitions.h"

PerfectHashFunction::PerfectHashFunction(Configuration config, ULLONG data_length, ULLONG *data) {
    ULLONG **bucket_data; //= new ULLONG*[_m]
    ULLONG *bucket_sizes; //= Array of ni's
    ULLONG max_bucket_size, max_mi;
    ULLONG *acyclicity_test_array;
    bool badTables, badFactor;
    short num_of_tries_tab = 0, num_of_tries_si;

    // RNG begin
    mt19937* rng;
    uniform_int_distribution<ULLONG>* dist_h_split_coeffs, dist_h_coeffs, dist_tables;
    // RNG end

    _configure(config, data_length);

    // RNG begin
    rng = new mt19937(config.seed);
    dist_h_split_coeffs = new uniform_int_distribution<ULLONG>(0, _h_split_mod_mask);
    // RNG end

    _h_split_coeffs = new ULLONG[_l+1];
    do {
        _createUhf(_h_split_coeffs, rng, dist_h_split_coeffs); //step 3
    } while(!_split(config, data_length, data, bucket_data, bucket_sizes, &max_bucket_size, &max_mi));
    //TODO assert(ceil(log(max_mi)) <= 64)

    // RNG begin
    dist_h_coeffs = new uniform_int_distribution<ULLONG>(0, _h_mod_mask);
    // RNG end

    _createGoodPairs(bucket_data, bucket_sizes, rng, dist_h_coeffs);

    // RNG begin
    dist_tables = new uniform_int_distribution<ULLONG>(0, (ULLONG)pow(2.0, _tab_width)-1);
    // RNG end

    acyclicity_test_array = new ULLONG[max_mi * 3];
    do {
        _createRandomTables(max_bucket_size, rng, dist_tables);
        badTables = false;
        num_of_tries_tab++;

        for(ULLONG i = 0; i < _m; i++) {
            num_of_tries_si = 0;
            do {
                _createRandomFactor(i, rng, dist_tables);
                num_of_tries_si++;
                _computeFij(i, acyclicity_test_array, bucket_sizes[i], bucket_data[i]);
                _computeGij(acyclicityTestArray);
                badFactor = _isCyclic(i, acyclicity_test_array);
            } while(badFactor && num_of_tries_si < config.num_of_tries_random_si);

            if(badFactor) {
                badTables = true;
                break;
            }
        }
    } while(badTables && num_of_tries_tab < config.num_of_tries_random_tab);

    if(badTables) {
        //TODO throw an exception!!!!
    }
}

void PerfectHashFunction::_configure(Configuration config, ULLONG data_length) {
    //TODO check this method!
    _k = config.k;
    _l = config.l;
    _m = (ULLONG)ceil(config.m_coeff * pow(data_length, config.m_exp));
    _h_split_mod_mask = (ULLONG)pow(2.0, _k + ceil(log2(_m)) + config.additional_bits_uhf)-1;
    _tab_rows = (ULLONG)ceil(config.tab_rows_coeff * pow(data_length, config.tab_rows_exp));
    _h_mod_mask = (ULLONG)pow(2.0, _k + ceil(log2(_tab_rows)) + config.additional_bits_uhf)-1;
}

void PerfectHashFunction::_createUhf(ULLONG* coeffs, mt19937* rng, uniform_int_distribution<ULLONG>* dist) {
    //TODO check this method!
    for(int i = 0; i < _l+1; i++) {
        coeffs[i] = (*dist)(*rng);
    }
}

ULLONG _evalUhf(ULLONG* coeff, ULLONG key) {
    ULLONG res = coeff[_l] & _h_split_mod_mask;
    // TODO
}

bool PerfectHashFunction::_split(Configuration config, data_length, ULLONG *data, ULLONG **bucket_data, ULLONG *bucket_sizes,
                                 ULLONG *max_bucket_size, ULLONG *max_mi) {
    //TODO check this method!
//    ULLONG **bucket_data; //= new ULLONG*[_m] // == B in concept
//    ULLONG *bucket_sizes; //= Array of ni's // <= diffs of C in concept
//    ULLONG max_bucket_size, max_mi;

    // initialization and instantiation
    ULLONG* splitted_data = new ULLONG[data_length]; // == B
    ULLONG* counter = new ULLONG[_m+1](); // initialized to zero because of (); == C[i]
    bucket_sizes = new ULLONG[_m]; // == n_i
    bucket_data = new ULLONG*[_m]; // bucket_data[i] == S[i] TODO or _m+1? (should not be needed...)
    ULLONG hv; // hash value by h_split
    ULLONG bucketOverflowSize = (ULLONG)floor(sqrt(data_length));

    // counting
    for(int i = 0; i < data_length; i++) {
        // TODO
    }
    // TODO


    // TODO step 5-7 + _tab_width

    return true;
}

void PerfectHashFunction::_createGoodPairs(ULLONG **bucket_data, ULLONG *bucket_sizes, mt19937* rng, uniform_int_distribution<ULLONG>* dist) {
    //TODO implement this method!
}

void PerfectHashFunction::_createRandomTables(ULLONG max_bucket_size, mt19937* rng, uniform_int_distribution<ULLONG>* dist) {
    //TODO implement this method!
}

void PerfectHashFunction::_createRandomFactor(ULLONG bucket_num, mt19937* rng, uniform_int_distribution<ULLONG>* dist) {
    //TODO implement this method!
}

void PerfectHashFunction::_computeFij(ULLONG bucket_num, ULLONG *acyclicity_test_array, ULLONG bucket_size,
                                      ULLONG *bucket) {
    //TODO implement this method!
}

void PerfectHashFunction::_computeGij(ULLONG *acyclicity_test_array) {
    //TODO implement this method!
}

bool PerfectHashFunction::_isCyclic(ULLONG bucket_num, ULLONG *acyclicity_test_array) {
    //TODO implement this method!
    return false;
}