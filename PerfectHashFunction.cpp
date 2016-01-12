//
// Created by philipp on 11.01.16.
//

#include "Definitions.h"

PerfectHashFunction::PerfectHashFunction(Configuration config, ULLONG data_length, ULLONG *data) {
    ULLONG **bucket_data; //= new ULLONG*[_m]
    ULLONG *bucket_sizes; //= Array of ni's
    ULLONG max_bucket_size, max_mi;
    ULLONG *acyclicity_test_array;
    bool badTables, badFactor;
    short num_of_tries_tab = 0, num_of_tries_si;

    _configure(config, data_length);
    do {
        _createUhf(_h_split_mod_mask, _h_split_coeffs); //step 3
    } while(!_split(data_length, data, bucket_data, bucket_sizes, &max_bucket_size, &max_mi));
    //TODO assert(ceil(log(max_mi)) <= 64)
    _createGoodPairs(bucket_data, bucket_sizes);
    acyclicity_test_array = new ULLONG[max_mi * 3];
    do {
        _createRandomTables(max_bucket_size);
        badTables = false;
        num_of_tries_tab++;

        for(ULLONG i = 0; i < _m; i++) {
            num_of_tries_si = 0;
            do {
                _createRandomFactor(i);
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

void PerfectHashFunction::_createUhf(ULLONG max_value, ULLONG* coeffs) {
    //TODO implement this method!
}

bool PerfectHashFunction::_split(ULLONG data_length, ULLONG *data, ULLONG **bucket_data, ULLONG *bucket_sizes,
                                 ULLONG *max_bucket_size, ULLONG *max_mi) {
    //TODO implement this method!
    return true;
}

void PerfectHashFunction::_createGoodPairs(ULLONG **bucket_data, ULLONG *bucket_sizes) {
    //TODO implement this method!
}

void PerfectHashFunction::_createRandomTables(ULLONG max_bucket_size) {
    //TODO implement this method!
}

void PerfectHashFunction::_createRandomFactor(ULLONG bucket_num) {
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