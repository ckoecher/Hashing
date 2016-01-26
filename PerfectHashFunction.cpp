//
// Created by philipp on 11.01.16.
//

#include "definitions.h"

PerfectHashFunction::PerfectHashFunction(Configuration config, ULLONG data_length, ULLONG *data) {
    ULLONG *splitted_data; // => B in (2)
    ULLONG **bucket_data; //= new ULLONG*[_m] => S[i] in (7) // TODO only needs normal delete[] because of splitted_data (!?)
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
    } while(!_split(config, data_length, data, splitted_data, bucket_data, bucket_sizes, &max_bucket_size, &max_mi));
    // TODO assert(ceil(log(max_mi)) <= 64)

    // RNG begin
    dist_h_coeffs = new uniform_int_distribution<ULLONG>(0, _h_mod_mask);
    // RNG end

    _createGoodPairs(bucket_data, bucket_sizes, max_bucket_size, rng, dist_h_coeffs);

    // RNG begin
    dist_tables = new uniform_int_distribution<ULLONG>(0, (ULLONG)pow(2.0, _tab_width)-1);
    // RNG end

    acyclicity_test_array = new ULLONG[max_mi * 3]();

    _random_table = new ULLONG[6*_tab_rows];
    _random_factor = new ULLONG[_m]; // TODO 3*_m for factors s_(i,j)
    do {
        _createRandomTables(rng, dist_tables);
        badTables = false;
        num_of_tries_tab++;

        for(ULLONG i = 0; i < _m; i++) {
            num_of_tries_si = 0;
            do {
                _createRandomFactor(i, rng, dist_tables);
                num_of_tries_si++;
                _computeGij(i, acyclicity_test_array, bucket_sizes[i], bucket_data[i]);
                badFactor = _isCyclic(i, acyclicity_test_array, bucket_sizes[i]);
            } while(badFactor && num_of_tries_si < config.num_of_tries_random_si);

            if(badFactor) {
                badTables = true;
                break;
            }
        }
    } while(badTables && num_of_tries_tab < config.num_of_tries_random_tab);

    if(badTables) {
        //TODO throw an exception!!!!
        // TODO deletes
    }
    // TODO deletes of short-term news
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

ULLONG _evalUhf(ULLONG key, ULLONG* coeff, ULLONG modMask, ULLONG modulus) {
    // TODO check this method!
    // & modMask == Mod (modMask+1) = Mod 2^l
    // modMask < 2^64
    // => x*y Mod 2^l == (x*y Mod 2^64) Mod 2^l
    // x*y Mod 2^64 == unsigned long long multiplication ("without" overflow)
    ULLONG keyMask = (ULLONG)pow(2.0, _k)-1;
    ULLONG res = coeff[_l] & modMask;
    for(int i = _l-1; i >= 0; i--) {
        res += (coeff[i] * (key & keyMask)) & modMask;
        res &= modMask;
        key >>= _k;
    }
    res >>= _k;
    res %= modulus;
    return res;
}

bool PerfectHashFunction::_split(Configuration config, ULLONG data_length, ULLONG* data, ULLONG* splitted_data, ULLONG **bucket_data, ULLONG *bucket_sizes,
                                 ULLONG *max_bucket_size, ULLONG *max_mi) {
    //TODO check this method!
//    ULLONG **bucket_data; //= new ULLONG*[_m] // == B in concept
//    ULLONG *bucket_sizes; //= Array of ni's // <= diffs of C in concept
//    ULLONG max_bucket_size, max_mi;

    // initialization and instantiation
    // ULLONG* splitted_data; // == B // external for better deallocation
    ULLONG* bucket_offsets; // == C[i]
    bucket_sizes = new ULLONG[_m](); // == n_i; initialized to zero because of ()
    ULLONG hv; // hash value by h_split
    ULLONG mi_1; // for computation of max_mi and _offset
    ULLONG bucketOverflowSize = (ULLONG)floor(sqrt(data_length)); // data_length >= 2^53 => maybe not exact result

    // counting
    for(int i = 0; i < data_length; i++) {
        hv = _evalUhf(data[i], _h_split_coeffs, _h_split_mod_mask, _m); // TODO save for usage in "sort data"?
        if(bucket_sizes[hv] >= bucketOverflowSize) {
            // bucket i overflow (too large)
            delete[] bucket_sizes;
            return false;
        } else {
            bucket_sizes[hv]++;
        }
    }
    // compute max. bucket size and offsets (for next segments) for sorted data (bucket_data segments)
    bucket_offsets = new ULLONG[_m+1](); // initialized to zero because of ()
    bucket_offsets[0] = bucket_sizes[0]; // TODO assert [0] exists?
    *max_bucket_size = bucket_sizes[0];
    for(int i = 1; i < _m; i++) {
        bucket_offsets[i] = bucket_offsets[i-1] + bucket_sizes[i];
        if(*max_bucket_size < bucket_sizes[i]) {
            *max_bucket_size = bucket_sizes[i];
        }
    }
    bucket_offsets[_m] = bucket_offsets[_m-1]; // TODO = data_length
    // sort data
    splitted_data = new ULLONG[data_length];
    for(int i = data_length-1; i >= 0; i--) {
        hv = _evalUhf(data[i], _h_split_coeffs, _h_split_mod_mask, _m);
        splitted_data[bucket_offsets[hv]-1] = data[i];
        bucket_offsets[hv]--;
    }
    bucket_data = new ULLONG*[_m]; // bucket_data[i] == S[i] TODO or _m+1? (should not be needed...)
    bucket_data[i] = splitted_data;
    for(int i = 1; i < _m; i++) {
        bucket_data[i] = bucket_data[i-1] + bucket_sizes[i-1];
    }
    // TODO _tab_width here or in _createRandomTables?
    _tab_width = (short)ceil(log2(*max_bucket_size))+(short)config.additional_bits_tab;
    // TODO assert(_tab_width <= 64);
    _offset = new ULLONG[_m+1];
    _offset[0] = 0;
    *max_mi = 0;
    for(int i = 1; i < _m+1; i++) {
        mi_1 = (ULLONG) ceil(config.mi_coeff * bucket_sizes[i-1]);
        if(*max_mi < mi_1) {
            *max_mi = mi_1;
        }
        _offset[i] = _offset[i-1] + mi_1;
    }

    // deallocate
    delete[] bucket_offsets;
    return true;
}

void PerfectHashFunction::_createGoodPairs(ULLONG **bucket_data, ULLONG *bucket_sizes, ULLONG max_bucket_size, mt19937* rng, uniform_int_distribution<ULLONG>* dist) {
    //TODO check this method!
    char* hTables = new char[(2*_tab_rows)/4+1]();
    ULLONG* hashValues = new ULLONG[2*max_bucket_size];
    ULLONG* h0coeffs, h1coeffs;
    bool goodPair;
    ULLONG x;

    _h_coeffs = new ULLONG[2*_m*(_l+1)];
    for(ULLONG pairI = 0; pairI < _m; pairI++) {
        // TODO Makro?
        h0coeffs = _h_coeffs + 2 * pairI * (_l + 1);
        h1coeffs = h0coeffs + _l + 1;
        goodPair = true;
        do {
            _createUhf(h0coeffs, rng, dist);
            _createUhf(h1coeffs, rng, dist);
            // counting
            for(ULLONG j = 0; j < bucket_sizes[i]; j++) {
                // compute hashvalues
                x = bucket_data[pairI][j];
                ARR(hashValues, bucket_sizes[i], 2, j, 0) = _evalUhf(x, h0coeffs, _h_mod_mask, _tab_width);
                ARR(hashValues, bucket_sizes[i], 2, j, 1) = _evalUhf(x, h1coeffs, _h_mod_mask, _tab_width);
                // count hashvalues
                if(GETBITPAIR(hTables, 0, ARR(hashValues, bucket_sizes[i], 2, j, 0)) < 2) {
                    INCBITPAIR(hTables, 0, ARR(hashValues, bucket_sizes[i], 2, j, 0));
                }
                if(GETBITPAIR(hTables, 1, ARR(hashValues, bucket_sizes[i], 2, j, 1)) < 2) {
                    INCBITPAIR(hTables, 1, ARR(hashValues, bucket_sizes[i], 2, j, 1));
                }
            }
            // evaluation
            for(ULLONG j = 0; j < bucket_sizes[i] && goodPair; j++) {
                // TODO &&
                if (GETBITPAIR(hTables, 0, ARR(hashValues, bucket_sizes[i], 2, j, 0)) == 2) {
                    if (GETBITPAIR(hTables, 1, ARR(hashValues, bucket_sizes[i], 2, j, 1)) == 2) {
                        // pair not good
                        goodPair = false;
                    }
                }
            }
            // clear hTables for counting
            // TODO without condition?
            if(pairI < _m-1 || !goodPair) {
                for(ULLONG j = 0; j < bucket_sizes[i]; j++) {
                    ZEROBITPAIRS(hTables, 0, ARR(hashValues, bucket_sizes[i], 2, j, 0));
                    ZEROBITPAIRS(hTables, 1, ARR(hashValues, bucket_sizes[i], 2, j, 1));
                }
            }
        } while(!goodPair);
    }
    delete[] hTables;
    delete[] hashValues;
}

void PerfectHashFunction::_createRandomTables(mt19937* rng, uniform_int_distribution<ULLONG>* dist) {
    //TODO check this method!
    for(int i = 0; i < 6*_tab_rows; i++) {
        _random_table[i] = (*dist)(*rng);
    }
}

void PerfectHashFunction::_createRandomFactor(ULLONG bucket_num, mt19937* rng, uniform_int_distribution<ULLONG>* dist) {
    //TODO check this method!
    _random_factor[bucket_num] = (*dist)(*rng);
}

void PerfectHashFunction::_computeGij(ULLONG bucket_num, ULLONG *acyclicity_test_array, ULLONG bucket_size,
                                      ULLONG *bucket) {
    //TODO check this method!
    ULLONG h0value, h1value, fi0, fi1, fi2;
    ULLONG *h0coeffs = _h_coeffs + 2 * bucket_num * (_l + 1);
    ULLONG *h1coeffs = h0coeffs + _l + 1;
    ULLONG x;
    ULLONG mi = _offset[bucket_num + 1] - _offset[bucket_num];
    for(ULLONG k = 0; k < bucket_size; k++) {
        x = bucket[k];
        h0value = _evalUhf(x, h0coeffs, _h_mod_mask, _tab_width);
        h1value = _evalUhf(x, h1coeffs, _h_mod_mask, _tab_width);

        //compute the values of the fij(x)
        fi0 = ((ARR(_random_table, _tab_rows, 6, h0value, 0) * _random_factor[bucket_num])
               ^ ARR(_random_table, _tab_rows, 6, h1value, 1)) % mi;
        fi1 = ((ARR(_random_table, _tab_rows, 6, h0value, 2) * _random_factor[bucket_num])
               ^ ARR(_random_table, _tab_rows, 6, h1value, 3)) % (mi - 1);
        fi2 = ((ARR(_random_table, _tab_rows, 6, h0value, 4) * _random_factor[bucket_num])
               ^ ARR(_random_table, _tab_rows, 6, h1value, 5)) % (mi - 2);

        //compute the values of the gij(x)
        if(fi1 >= fi0) {
            fi1++;
        }
        if(fi2 >= fi0) {
            if(fi2 >= fi1) {
                fi2 += 2;
            } else {
                fi2++;
            }
        } else if(fi2 >= fi1) {
            fi2++;
        }

        //save the values in the array
        ARR(acyclicity_test_array, bucket_size, 3, k, 0) = fi0;
        ARR(acyclicity_test_array, bucket_size, 3, k, 1) = fi1;
        ARR(acyclicity_test_array, bucket_size, 3, k, 2) = fi2;
    }
}

bool PerfectHashFunction::_isCyclic(ULLONG bucket_num, ULLONG *acyclicity_test_array, ULLONG bucket_size) {
    ULLONG max_length = 2*log2(bucket_size);
    ULLONG mi = _offset[bucket_num + 1] - _offset[bucket_num];
    ULLONG *edgesOf = new ULLONG[max_length * mi]();
    ULLONG *cEdgesOf = new ULLONG[mi]();
    ULLONG gValue;
    bool *removed = new bool[bucket_size](); //TODO shouldn't we use char, int or something similar here?
    std::queue<ULLONG> queue;

    //construct a list of adjacent edges of each node
    for(ULLONG j = 0; j < bucket_size; j++) {
        for(int k = 0; k < 3; k++) {
            gValue = ARR(acyclicity_test_array, bucket_size, 3, j, k);
            ARR(edgesOf, mi, max_length, gValue, cEdgesOf[gValue]) = j;
            cEdgesOf[gValue]++;
            if(cEdgesOf[gValue] == max_length) { //an overflow here
                delete[] edgesOf;
                delete[] cEdgesOf;
                delete[] removed;
                delete[] queue;
                return true; //TODO is this right?
            }
        }
    }

    //now check for acyclicity
    for(ULLONG j = 0; j < mi; j++) {
        if(cEdgesOf[j] == 1 && !removed[ARR(edgesOf, mi, max_length, j, 0)]) {
            peelOf(0, j, queue, edgesOf, cEdgesOf, max_length, mi);
        }
    }

    if(queue.size() != bucket_size) {
        delete[] edgesOf;
        delete[] cEdgesOf;
        delete[] removed;
        delete[] queue;
        return true;
    }

    //now assign the values
    //TODO implement this method!

    delete[] edgesOf;
    delete[] cEdgesOf;
    delete[] removed;
    delete[] queue;
    return false;
}

void peelOf(ULLONG edge_index, ULLONG vertex_index, std::queue<ULLONG> queue, ULLONG *edgesOf, ULLONG *cEdgesOf,
            ULLONG max_length, ULLONG mi) {
    //TODO implement this method!
}