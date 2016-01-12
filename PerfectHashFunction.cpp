//
// Created by philipp on 11.01.16.
//

#include "Definitions.h"

PerfectHashFunction::PerfectHashFunction(Configuration config, ULLONG data_length, ULLONG *data) {
    ULLONG **bucket_data; //= new ULLONG*[_m]
    ULLONG *bucket_sizes; //= Array of ni's
    ULLONG max_bucket_size, max_mi;
    ULLONG *acyclicityTestArray;
    bool badTables, badFactor;
    short num_of_tries_tab = 0, num_of_tries_si;

    _configure(config, data_length);
    do {
        _createUhf(_h_split_mod_mask, _h_split_coeffs); //step 3
    } while(!_split(data_length, data, bucket_data, bucket_sizes, &max_bucket_size, &max_mi));
    //TODO assert(ceil(log(max_mi)) <= 64)
    _createGoodPairs(bucket_data, bucket_sizes);
    _acyclicityTestArray = new ULLONG[max_mi * 3];
    do {
        _createRandomTables(max_bucket_size);
        badTables = false;
        num_of_tries_tab++;

        for(ULLONG i = 0; i < _m; i++) {
            num_of_tries_si = 0;
            do {
                _createRandomFactor(i);
                num_of_tries_si++;
                _computeFij(i, acyclicityTestArray, bucket_sizes[i], bucket_data[i]);
                _computeGij(acyclicityTestArray);
                badFactor = _isCyclic(i, acyclicityTestArray);
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