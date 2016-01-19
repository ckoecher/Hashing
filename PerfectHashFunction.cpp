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
    ULLONG maxRandomNumber;
    // RNG end

    _configure(config, data_length);

    // RNG begin
    rng = new mt19937(config.seed);
    maxRandomNumber = (ULLONG)pow(2.0, _k + ceil(log2(_m)) + additional_bits_uhf)-1;
    dist_h_split_coeffs = new uniform_int_distribution<ULLONG>(0,maxRandomNumber);
    // RNG end

    do {
        _createUhf(_h_split_mod_mask, _h_split_coeffs, rng, dist_h_split_coeffs); //step 3
    } while(!_split(data_length, data, bucket_data, bucket_sizes, &max_bucket_size, &max_mi));
    //TODO assert(ceil(log(max_mi)) <= 64)

    // RNG begin
    maxRandomNumber = (ULLONG)pow(2.0, _k + ceil(log2(_tab_rows)) + config.additional_bits_uhf)-1;
    dist_h_coeffs = new uniform_int_distribution<ULLONG>(0,maxRandomNumber);
    // RNG end

    _createGoodPairs(bucket_data, bucket_sizes, rng, dist_h_coeffs);

    // RNG begin
    maxRandomNumber = (ULLONG)pow(2.0, ceil(log2(max_bucket_size)) + config.additional_bits_tab)-1;
    dist_tables = new uniform_int_distribution<ULLONG>(0,maxRandomNumber);
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
    //TODO implement this method!
}

void PerfectHashFunction::_createUhf(ULLONG max_value, ULLONG* coeffs, mt19937* rng, uniform_int_distribution<ULLONG>* dist) {
    //TODO implement this method!
}

bool PerfectHashFunction::_split(ULLONG data_length, ULLONG *data, ULLONG **bucket_data, ULLONG *bucket_sizes,
                                 ULLONG *max_bucket_size, ULLONG *max_mi) {
    //TODO implement this method!
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