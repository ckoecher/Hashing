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

    configure(config, data_length);
    do {
        createUhf(_h_split_mod_mask, _h_split_coeffs); //step 3
    } while(!split(data_length, data, bucket_data, bucket_sizes, &max_bucket_size, &max_mi));
    //TODO assert(ceil(log(max_mi)) <= 64)
    createGoodPairs(bucket_data, bucket_sizes);
    acyclicityTestArray = new ULLONG[max_mi * 3];
    do {
        createRandomTables(max_bucket_size);
        badTables = false;
        num_of_tries_tab++;

        for(ULLONG i = 0; i < _m; i++) {
            num_of_tries_si = 0;
            do {
                createRandomFactor(i);
                num_of_tries_si++;
                computeFij(i, acyclicityTestArray, bucket_sizes[i], bucket_data[i]);
                computeGij(acyclicityTestArray);
                badFactor = isCyclic(i, acyclicityTestArray);
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