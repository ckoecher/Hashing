//
// Created by philipp on 11.01.16.
//

#include "definitions.h"
#include "PerfectHashFunction.h"

PerfectHashFunction::PerfectHashFunction(Configuration &config, InputData *data, Statistics &stats) {
    // TODO check this implementation

    // Stats
    stats.creation_start = clock();
    stats.num_of_keys = data->getLength();
    stats.setup_start = clock();
    // Stats end

    InputData *bucket_data = new InputData();
    ULLONG *bucket_offsets = nullptr;
    ULLONG max_bucket_size, max_mi;
    ULLONG *acyclicity_test_array = nullptr;
    bool badTables, badFactor;
    bool split_successful;
    unsigned short num_of_tries_split = 0;
    unsigned short num_of_tries_tab = 0, num_of_tries_si;

    // Debug
    if(config.debug_mode) {
        cout << "### Perfect Hash Function Creation ###" << endl;
        cout << data->getLength() << " key(s)." << endl;
//        cout << "data[i]:";
//        for(ULLONG i = 0; i < data->getLength(); i++) {
//            cout << " " << data->getValue(i);
//        }
//        cout << endl;
        cout << "# Setup..." << endl;
    }
    // Debug end

    // RNG begin
    mt19937* rng = nullptr;
    uniform_int_distribution<ULLONG> *dist_h_split_coeffs = nullptr, *dist_h_coeffs = nullptr, *dist_tables = nullptr;
    // RNG end

    _configure(config, data->getLength());

    // RNG begin
    rng = new mt19937(config.seed);
    dist_h_split_coeffs = new uniform_int_distribution<ULLONG>(0, _h_split_mod_mask);
    // RNG end

    // Debug
    if(_debug_mode) {
        cout << "Setup successful." << endl;
        cout << "# Split data into buckets:" << endl;
    }
    // Debug end

    // Stats
    stats.setup_end = clock();
    stats.setup_time = stats.setup_end - stats.setup_start;
    stats.setup_succuess = true;
    stats.split_start = clock();
    // Stats end

    _h_split_coeffs = new ULLONG[_l+1];
    do {
        num_of_tries_split++;
        // Debug
        if(_debug_mode) {
            cout << "Split no. " << num_of_tries_split << endl;
        }
        // Debug end
        _createUhf(_h_split_coeffs, rng, dist_h_split_coeffs); //step 3
        split_successful = _split(config, data, bucket_data, bucket_offsets, &max_bucket_size, &max_mi, stats);
    } while(!split_successful && num_of_tries_split < config.num_of_tries_split);
    // TODO assert(ceil(log(max_mi)) <= 64)

    if(!split_successful) {
        // TODO
        // Debug
        if(_debug_mode) {
            cout << "Could not split data within " << num_of_tries_split << " tries." << endl;
        }
        // Debug end
        throw 0;
    }

    // Debug
    if(_debug_mode) {
        cout << "Could split data within " << num_of_tries_split << " tries." << endl;
        cout << "Created " << _m << " buckets with <= " << max_bucket_size << " elements." << endl;
        cout << "Bucket sizes:" << endl;
        for(ULLONG i = 0; i < _m; i++) {
            cout << " " << bucket_offsets[i+1]-bucket_offsets[i];
        }
        cout << endl;
    }
    // Debug end

    // Stats
    stats.split_end = clock();
    stats.split_time = stats.split_end - stats.split_start;
    stats.split_success = true;
    stats.split_tries = num_of_tries_split;
    stats.num_of_buckets = _m;
    stats.max_bucket_size = max_bucket_size;
    //stats.min_bucket_size computed in _split(...)
    stats.avg_bucket_size = (long double)stats.num_of_keys / (long double)stats.num_of_buckets;
    stats.goodpairs_start = clock();
    // Stats end

//    // Debug
//    cout << "### after split ###" << endl;
//    cout << "bucket_sizes[i]:" << endl;
//    for(ULLONG i = 0; i < _m; i++) {
//        cout << " " << bucket_offsets[i+1]-bucket_offsets[i];
//    }
//    // TODO bucket_sizes[i] = bucket_offsets[i+1]-bucket_offsets[i]
//    cout << endl;
//    cout << "bucket_offsets[i]:" << endl;
//    for(ULLONG i = 0; i < _m+1; i++) {
//        cout << " " << bucket_offsets[i];
//    }
//    cout << endl;
//    cout << "bucket_data[i]:" << endl;
//    for(ULLONG i = 0; i < _m; i++) {
//        cout << "\tbucket " << i << ":";
//        for(ULLONG j = 0; j < bucket_offsets[i+1]-bucket_offsets[i]; j++) {
//            //cout << " " << bucket_data[i][j];
//            cout << " " << bucket_data->getValue(bucket_offsets[i]+j);
//        }
//        cout << endl;
//    }
//    // Debug end

    // RNG begin
    dist_h_coeffs = new uniform_int_distribution<ULLONG>(0, _h_mod_mask);
    // RNG end

    // Debug
    if(_debug_mode) {
        cout << "# Create good pairs of universal hash functions:" << endl;
    }
    // Debug end

    _createGoodPairs(config, bucket_data, bucket_offsets, max_bucket_size, rng, dist_h_coeffs, stats);

//    // Debug
//    cout << "### \"testing\" good pairs ###" << endl;
//    cout << "bucket i: (h_0^i(x), h_1^i(x)) (h_0^i(y), h_1^i(y)) ..." << endl;
//    for(ULLONG i = 0; i < _m; i++) {
//        cout << "bucket " << i << ":";
//        ULLONG *h0coeffs = _h_coeffs + ((i * (_l + 1)) << 1);
//        ULLONG *h1coeffs = h0coeffs + _l + 1;
//        for(int j = 0; j < bucket_offsets[i+1] - bucket_offsets[i]; j++) {
//            cout << " (" << _evalUhf(bucket_data->getValue(bucket_offsets[i]+j), h0coeffs, _h_mod_mask, _tab_rows)
//                 << ", " << _evalUhf(bucket_data->getValue(bucket_offsets[i]+j), h1coeffs, _h_mod_mask, _tab_rows)
//                 << ")";
//        }
//        cout << endl;
//    }
//    // Debug end

    // Debug
    if(_debug_mode) {
        cout << "Created good pairs of universal hash functions for all " << _m << " buckets." << endl;
    }
    // Debug end

    // Stats
    stats.goodpairs_end = clock();
    stats.goodpairs_time = stats.goodpairs_end - stats.goodpairs_start;
    stats.goodpairs_success = true;
    stats.buckets_start = clock();
    // Stats end

    // RNG begin
    dist_tables = new uniform_int_distribution<ULLONG>(0, (ULLONG)pow(2.0l, _tab_width)-1);
    // RNG end

    // Debug
    if(_debug_mode) {
        cout << "# Create perfect hash functions for all buckets:" << endl;
    }
    // Debug end

    acyclicity_test_array = new ULLONG[max_mi * 3]();

    _random_table = new ULLONG[6*_tab_rows];
    _random_factor = new ULLONG[_m]; // TODO 3*_m for factors s_(i,j)
    _g = new unsigned char[(_offset[_m] >> 2) + 1]();

    // TODO debug following code

    do {

        num_of_tries_tab++;
        // Stats
        stats.random_tab_tries++;
        // Stats end

        // Debug
        if(_debug_mode) {
            cout << "Create new random tables (no. " << num_of_tries_tab << ")." << endl;
        }
        // Debug end

        _createRandomTables(rng, dist_tables);
        badTables = false;

        for(ULLONG i = 0; i < _m; i++) {
//            // Debug
//            cout << "# bucket " << i << " with " << bucket_offsets[i+1]-bucket_offsets[i] << " elements #" << endl;
//            for(ULLONG jj = 0; jj < bucket_offsets[i+1]-bucket_offsets[i]; jj++) {
//                cout << " " << bucket_data->getValue(bucket_offsets[i]+jj) << endl;
//            }
//            // Debug end
            // Debug
            if(_debug_mode) {
                cout << "Bucket " << i << ":" << endl;
            }
            // Debug end
            num_of_tries_si = 0;
            do {
                num_of_tries_si++;
                // Stats
                stats.random_si_total_tries++;
                // Stats end
                // Debug
                if(_debug_mode) {
                    cout << " Create new random factor (no. " << num_of_tries_si << ") for bucket " << i << "." << endl;
                }
                // Debug end
                _createRandomFactor(i, rng, dist_tables);
                // Debug
                if(_debug_mode) {
                    cout << " Create 3-graph." << endl;
                }
                // Debug end
                _computeGij(i, acyclicity_test_array, bucket_data, bucket_offsets[i], bucket_offsets[i+1]-bucket_offsets[i]);
                // Debug
                if(_debug_mode) {
                    cout << " Test 3-graph for acyclicity." << endl;
                }
                // Debug end
                badFactor = _isCyclic(i, acyclicity_test_array, bucket_offsets[i+1]-bucket_offsets[i]);
            } while(badFactor && num_of_tries_si < config.num_of_tries_random_si);

            if(badFactor) {
                // Debug
                if(_debug_mode) {
                    cout << " Could not create acyclic 3-graph for bucket " << i << " within " << num_of_tries_si << " tries." << endl;
                }
                // Debug end
                badTables = true;
                break;
            }
        }
    } while(badTables && num_of_tries_tab < config.num_of_tries_random_tab);

    // Stats
    stats.buckets_end = clock();
    stats.buckets_time = stats.buckets_end - stats.buckets_start;
    stats.buckets_success = !badTables;
    // Stats end

    // delete temporary data
    // delete nullptr has no effect
    //delete[] splitted_data;
    bucket_data->close();
    delete bucket_data;
    //delete[] bucket_sizes;
    delete[] acyclicity_test_array;
    delete rng;
    delete dist_h_split_coeffs;
    delete dist_h_coeffs;
    delete dist_tables;

    if(badTables) {
        // Debug
        if(_debug_mode) {
            cout << "Could not create random tables within " << num_of_tries_tab << " tries." << endl;
        }
        // Debug end
        // construction not successful
        _clear();
        //TODO throw an exception!!!!
        throw 0;
    }

    // Stats
    stats.creation_end = clock();
    stats.creation_time = stats.creation_end - stats.creation_start;
    stats.creation_success = true;
    stats.range_of_phf = getRange();
    stats.size_in_bytes = getSizeInBytes();
    // Stats end

    // Debug
    if(_debug_mode) {
        cout << "### Perfect Hash Function Creation Successful ###" << endl;
    }
    // Debug end
}

void PerfectHashFunction::_configure(Configuration &config, ULLONG data_length) {
    _k = config.k;
    _l = config.l;
    _m = (ULLONG)ceil(config.m_coeff * pow(data_length, config.m_exp));
    _h_split_mod_mask = (ULLONG)pow(2.0l, _k + ceil(log2(_m)) + config.additional_bits_uhf)-1;
    _tab_rows = (ULLONG)ceil(config.tab_rows_coeff * pow(data_length, config.tab_rows_exp));
    _h_mod_mask = (ULLONG)pow(2.0l, _k + ceil(log2(_tab_rows)) + config.additional_bits_uhf)-1;
    _debug_mode = config.debug_mode;
}

void PerfectHashFunction::_createUhf(ULLONG* coeffs, mt19937* rng, uniform_int_distribution<ULLONG>* dist) {
    for(int i = 0; i < _l+1; i++) {
        coeffs[i] = (*dist)(*rng);
    }
}

ULLONG PerfectHashFunction::_evalUhf(ULLONG key, ULLONG* coeff, ULLONG modMask, ULLONG modulus) {
    // & modMask == Mod (modMask+1) = Mod 2^l
    // modMask < 2^64
    // => x*y Mod 2^l == (x*y Mod 2^64) Mod 2^l
    // x*y Mod 2^64 == unsigned long long multiplication ("without" overflow)
    ULLONG keyMask = (ULLONG)pow(2.0l, _k)-1;
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

bool PerfectHashFunction::_split(Configuration &config, InputData *data, InputData *bucket_data,
                                 ULLONG *& bucket_offsets, ULLONG *max_bucket_size, ULLONG *max_mi, Statistics &stats) {
    // initialization and instantiation
    ULLONG data_length = data->getLength(); // = n
    ULLONG* bucket_sizes = new ULLONG[_m](); // == n_i; initialized to zero because of ()
    ULLONG hv; // hash value by h_split
    ULLONG mi_1; // for computation of max_mi and _offset
    ULLONG bucketOverflowSize = (ULLONG)floor(sqrt(data_length)); // data_length >= 2^53 => maybe not exact result
    ULLONG value; // temporary value from data stream

//    // Debug
//    if(_debug_mode) {
//        cout << "bucket_sizes[i]:";
//        for(ULLONG i = 0; i < _m; i++) {
//            cout << " " << bucket_sizes[i];
//        }
//        cout << endl;
//    }
//    // Debug end

//    // Debug
//    cout << "buckets:";
//    // Debug end

    // counting
    for(ULLONG i = 0; i < data_length; i++) {
        hv = _evalUhf(data->getValue(i), _h_split_coeffs, _h_split_mod_mask, _m); // TODO save for usage in "sort data"?
//        // Debug
//        cout << " " << hv;
//        // Debug end
        if(bucket_sizes[hv] >= bucketOverflowSize) {
            // bucket i overflow (too large)

//            // Debug
//            cout << "_split not possible with _h_split_coeffs";
//            for(ULLONG jj = 0; jj < _l+1; jj++) {
//                cout << " " << _h_split_coeffs[jj];
//            }
//            cout << endl;
//            // Debug end

            delete[] bucket_sizes;
            return false;
        } else {
            bucket_sizes[hv]++;
        }
    }
//    // Debug
//    cout << endl;
//    // Debug end

//    // Debug
//    if(_debug_mode) {
//        cout << "bucket_sizes[i]:";
//        for(ULLONG i = 0; i < _m; i++) {
//            cout << " " << bucket_sizes[i];
//        }
//        cout << endl;
//    }
//    // Debug end

    // compute max. bucket size and offsets (for next segments) for sorted data (bucket_data segments)
    bucket_offsets = new ULLONG[_m+1](); // initialized to zero because of ()
    bucket_offsets[0] = bucket_sizes[0]; // TODO assert [0] exists?
    *max_bucket_size = bucket_sizes[0];
    // Stats
    stats.min_bucket_size = bucket_sizes[0];
    // Stats end
    for(int i = 1; i < _m; i++) {
        bucket_offsets[i] = bucket_offsets[i-1] + bucket_sizes[i];
        if(*max_bucket_size < bucket_sizes[i]) {
            *max_bucket_size = bucket_sizes[i];
        }
        // Stats
        if(stats.min_bucket_size > bucket_sizes[i]) {
            stats.min_bucket_size = bucket_sizes[i];
        }
        // Stats end
    }
    bucket_offsets[_m] = bucket_offsets[_m-1]; // TODO = data_length

/*    // Debug
    if(_debug_mode) {
        cout << "bucket_offsets[i]:";
        for(ULLONG i = 0; i < _m+1; i++) {
            cout << " " << bucket_offsets[i];
        }
        cout << endl;
    }
    // Debug end*/

    // sort data
    //splitted_data = new ULLONG[data_length];
    for(ULLONG i = data_length-1; i < data_length; i--) { // condition: not i >= 0 because of UNSIGNED
        value = data->getValue(i);
        hv = _evalUhf(value, _h_split_coeffs, _h_split_mod_mask, _m);
        //splitted_data[bucket_offsets[hv]-1] = value;
        bucket_data->setValue(value, bucket_offsets[hv]-1);
        bucket_offsets[hv]--;
    }
    //bucket_data = new ULLONG*[_m]; // bucket_data[i] == S[i] TODO or _m+1? (should not be needed...)
    //bucket_data[0] = splitted_data;
    //for(ULLONG i = 1; i < _m; i++) {
    //    bucket_data[i] = bucket_data[i-1] + bucket_sizes[i-1];
    //}
    // TODO _tab_width here or in _createRandomTables?
    _tab_width = (unsigned short)ceil(log2(*max_bucket_size))+config.additional_bits_tab;
    // TODO assert(_tab_width <= 64);
    _offset = new ULLONG[_m+1];
    _offset[0] = 0;
    *max_mi = 0;
    for(ULLONG i = 1; i < _m+1; i++) {
        mi_1 = (ULLONG) ceil(config.mi_coeff * bucket_sizes[i-1]);
        // TODO maybe not necessary (to eliminate problem with %(mi-2))
        if(mi_1 > 0 && mi_1 < 3) {
            mi_1 = 3;
        }
        if(*max_mi < mi_1) {
            *max_mi = mi_1;
        }
        _offset[i] = _offset[i-1] + mi_1;
    }

//    // Debug
//    if(_debug_mode) {
//        cout << "_offset[i]:";
//        for(ULLONG i = 0; i < _m+1; i++) {
//            cout << " " << _offset[i];
//        }
//        cout << endl;
//    }
//    // Debug end

    // deallocate
    delete[] bucket_sizes;
    return true;
}

void PerfectHashFunction::_createGoodPairs(Configuration &config, InputData *bucket_data, ULLONG *bucket_offsets, ULLONG max_bucket_size, mt19937* rng, uniform_int_distribution<ULLONG>* dist, Statistics &stats) {
    // TODO debug following code (if necessary)
    unsigned char *hTables = new unsigned char[(_tab_rows>>1)+1](); // (2*_tab_rows)/4+1
    ULLONG *hashValues = new ULLONG[max_bucket_size<<1]; // 2*max_bucket_size
    ULLONG *h0coeffs = nullptr, *h1coeffs = nullptr;
    bool goodPair;
    ULLONG x, bucket_size;
    unsigned short num_of_tries_goodpairs;

    _h_coeffs = new ULLONG[(_m*(_l+1))<<1]; // 2*_m*(_l+1)
    for(ULLONG pairI = 0; pairI < _m; pairI++) {

        num_of_tries_goodpairs = 0;
//        // Debug
//        cout << "## create good pair for bucket " << pairI << " ##" << endl;
//        // Debug end

        // TODO Makro?
        h0coeffs = _h_coeffs + ((pairI * (_l + 1))<<1); // _h_coeffs + 2 * pairI * (_l + 1)
        h1coeffs = h0coeffs + _l + 1;
        // goodPair = true;
        bucket_size = bucket_offsets[pairI+1]-bucket_offsets[pairI];

        // Debug
        if(_debug_mode) {
            cout << "Bucket " << pairI << ":" << endl;
            cout << " Need to assign " << bucket_size << " elements to " << _tab_rows << " different values." << endl;
        }
        // Debug end

        do {
            goodPair = true;
            num_of_tries_goodpairs++;
            // Stats
            stats.goodpairs_total_tries++;
            // Stats end

            _createUhf(h0coeffs, rng, dist);
            _createUhf(h1coeffs, rng, dist);
            // counting
            //for(ULLONG j = 0; j < bucket_sizes[pairI]; j++) {
            for(ULLONG j = 0; j < bucket_size; j++) {
                // compute hashvalues
                //x = bucket_data[pairI][j];
                x = bucket_data->getValue(bucket_offsets[pairI]+j);
                ARR(hashValues, bucket_size, 2, j, 0) = _evalUhf(x, h0coeffs, _h_mod_mask, _tab_rows);
                ARR(hashValues, bucket_size, 2, j, 1) = _evalUhf(x, h1coeffs, _h_mod_mask, _tab_rows);
                // count hashvalues
                if(GETBITPAIR(hTables, 0, ARR(hashValues, bucket_size, 2, j, 0)) < 2) {
                    INCBITPAIR(hTables, 0, ARR(hashValues, bucket_size, 2, j, 0));
                }
                if(GETBITPAIR(hTables, 1, ARR(hashValues, bucket_size, 2, j, 1)) < 2) {
                    INCBITPAIR(hTables, 1, ARR(hashValues, bucket_size, 2, j, 1));
                }
            }
            // evaluation
            for(ULLONG j = 0; j < bucket_size && goodPair; j++) {
                // TODO &&
                if (GETBITPAIR(hTables, 0, ARR(hashValues, bucket_size, 2, j, 0)) == 2) {
                    if (GETBITPAIR(hTables, 1, ARR(hashValues, bucket_size, 2, j, 1)) == 2) {
                        // pair not good
                        goodPair = false;
                    }
                }
            }
            // clear hTables for counting
            // TODO without condition?
            if(pairI < _m-1 || !goodPair) {
                for(ULLONG j = 0; j < bucket_size; j++) {
                    ZEROBITPAIRS(hTables, 0, ARR(hashValues, bucket_size, 2, j, 0));
                    ZEROBITPAIRS(hTables, 1, ARR(hashValues, bucket_size, 2, j, 1));
                }
            }

//            // Debug
//            if(!goodPair) {
//                cout << "rerun" << endl;
//            }
//            // Debug end
        } while(!goodPair && num_of_tries_goodpairs < config.num_of_tries_goodpairs);

        if(!goodPair) {
            // Debug
            if(_debug_mode) {
                cout << " Could not create good pair of universal hash functions for bucket " << pairI
                     << " within " << num_of_tries_goodpairs << " tries." << endl;
            }
            // Debug end
            // TODO
            throw 0;
        }

        // Debug
        if(_debug_mode) {
            cout << " Good pair " << pairI << " created within " << num_of_tries_goodpairs << " tries." << endl;
        }
        // Debug end
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

void PerfectHashFunction::_computeGij(ULLONG bucket_num, ULLONG *acyclicity_test_array,
                                      InputData *bucket_data, ULLONG bucket_offset, ULLONG bucket_size) {
    //TODO check this method!
    ULLONG h0value, h1value, fi0, fi1, fi2;
    ULLONG *h0coeffs = _h_coeffs + ((bucket_num * (_l + 1))<<1); // _h_coeffs + 2 * bucket_num * (_l + 1)
    ULLONG *h1coeffs = h0coeffs + _l + 1;
    ULLONG x;
    ULLONG mi = _offset[bucket_num + 1] - _offset[bucket_num];
    for(ULLONG k = 0; k < bucket_size; k++) {
        //x = bucket[k];
        x = bucket_data->getValue(bucket_offset+k);
        h0value = _evalUhf(x, h0coeffs, _h_mod_mask, _tab_rows);
        h1value = _evalUhf(x, h1coeffs, _h_mod_mask, _tab_rows);

        //compute the values of the fij(x)
        fi0 = ((ARR(_random_table, _tab_rows, 6, h0value, 0) * _random_factor[bucket_num])
               ^ (ARR(_random_table, _tab_rows, 6, h1value, 1))) % mi;
        fi1 = ((ARR(_random_table, _tab_rows, 6, h0value, 2) * _random_factor[bucket_num])
               ^ (ARR(_random_table, _tab_rows, 6, h1value, 3))) % (mi - 1);
        fi2 = ((ARR(_random_table, _tab_rows, 6, h0value, 4) * _random_factor[bucket_num])
               ^ (ARR(_random_table, _tab_rows, 6, h1value, 5))) % (mi - 2);

        //compute the values of the gij(x)
        if(fi1 >= fi0) {
            fi1++;
        }
//        if(fi2 >= fi0) {
//            if(fi2 >= fi1) {
//                fi2 += 2;
//            } else {
//                fi2++;
//            }
//        } else if(fi2 >= fi1) {
//            fi2++;
//        }
        // first version not equivalent / not correct
        if(fi2 >= fi0) {
            fi2++;
            if(fi2 >= fi1) {
                fi2++;
            }
        } else if(fi2 >= fi1) {
            fi2++;
            if(fi2 >= fi0) {
                fi2++;
            }
        }

        //save the values in the array
        ARR(acyclicity_test_array, bucket_size, 3, k, 0) = fi0;
        ARR(acyclicity_test_array, bucket_size, 3, k, 1) = fi1;
        ARR(acyclicity_test_array, bucket_size, 3, k, 2) = fi2;
    }
}

bool PerfectHashFunction::_isCyclic(ULLONG bucket_num, ULLONG *acyclicity_test_array, ULLONG bucket_size) {
    //TODO check this method!
    // TODO new for queue, removed and visited after first loop?
    ULLONG max_length = 2*ceil(log2(bucket_size))+1; // TODO <<1 ??? (after cast?)
    ULLONG mi = _offset[bucket_num + 1] - _offset[bucket_num];
    ULLONG *edgesOf = new ULLONG[max_length * mi]();
    ULLONG *cEdgesOf = new ULLONG[mi]();
    ULLONG gValue, edge_index, u;
    ULLONG *queue = new ULLONG[bucket_size];
    ULLONG next_queue_index = 0;
    unsigned char *removed = new unsigned char[(bucket_size >> 3) + 1]();
    unsigned char *visited = new unsigned char[(mi >> 3) + 1]();
    int c, sum;

//    // Debug
//    cout << "Check new 3-graph with " << bucket_size << " edges and "
//         << _offset[bucket_num+1]-_offset[bucket_num] << " nodes" << endl;
//    cout << "list of edges:" << endl;
//    for(ULLONG ii = 0; ii < bucket_size; ii++) {
//        cout << "\t";
//        for(ULLONG jj = 0; jj < 3; jj++) {
//            cout << ARR(acyclicity_test_array, bucket_offsets[i+1]-bucket_offsets[i], 3, ii, jj) << " ";
//        }
//        cout << endl;
//    }
//    // Debug end

    //construct a list of adjacent edges of each node
    for(ULLONG j = 0; j < bucket_size; j++) {
        for(int k = 0; k < 3; k++) {
            gValue = ARR(acyclicity_test_array, bucket_size, 3, j, k);
            if(cEdgesOf[gValue] >= max_length) { //an overflow here
//                // Debug
//                cout << "node " << gValue << " has more than " << max_length << " edges" << endl;
//                cout << "creation not successful" << endl;
//                // Debug end
                delete[] edgesOf;
                delete[] cEdgesOf;
                delete[] queue;
                delete[] removed;
                delete[] visited;
                return true; //TODO is this right?
            }
            ARR(edgesOf, mi, max_length, gValue, cEdgesOf[gValue]) = j;
            cEdgesOf[gValue]++;
        }
    }

//    // Debug
//    cout << "list of adjacent edges:" << endl;
//    for(ULLONG nodeindex = 0; nodeindex < mi; nodeindex++) {
//        cout << "\tnode " << nodeindex << " with " << cEdgesOf[nodeindex] << " edges:";
//        for(int edgeindex = 0; edgeindex < cEdgesOf[nodeindex]; edgeindex++) {
//            cout << " " << ARR(edgesOf, mi, max_length, nodeindex, edgeindex);
//        }
//        cout << endl;
//    }
//    // Debug end

    //now check for acyclicity
    for(ULLONG j = 0; j < mi; j++) {
        if(cEdgesOf[j] == 1) {
            edge_index = ARR(edgesOf, mi, max_length, j, 0);
//            // Debug
//            cout << "edge " << edge_index << " already removed? " << GETBIT(removed, edge_index) << endl;
//            cout << "repeat " << GETBIT(removed, edge_index) << endl;
//            // Debug end
            if(!GETBIT(removed, edge_index)) {
//                // Debug
//                cout << "remove edge " << edge_index << endl;
//                // Debug end
                _peelOf(edge_index, j, acyclicity_test_array, bucket_size, queue, next_queue_index, edgesOf,
                       cEdgesOf, max_length, mi, removed);
            }
        }
    }

    if(next_queue_index != bucket_size) {
//        // Debug
//        if(_debug_mode) {
//            cout << "bucket " << bucket_num << ": not acyclic" << endl;
//        }
//        // Debug end
        delete[] edgesOf;
        delete[] cEdgesOf;
        delete[] queue;
        delete[] removed;
        delete[] visited;
        return true;
    }

//    // Debug
//    if(_debug_mode) {
//        cout << "bucket " << bucket_num << ": acyclic" << endl;
//    }
//    // Debug end
//    // Debug
//    cout << "list of adjacent edges:" << endl;
//    for(ULLONG nodeindex = 0; nodeindex < mi; nodeindex++) {
//        cout << "\tnode " << nodeindex << " with " << cEdgesOf[nodeindex] << " edges:";
//        for(int edgeindex = 0; edgeindex < cEdgesOf[nodeindex]; edgeindex++) {
//            cout << " " << ARR(edgesOf, mi, max_length, nodeindex, edgeindex);
//        }
//        cout << endl;
//    }
//    // Debug end

    //now assign the values
    for(ULLONG u = _offset[bucket_num]; u < _offset[bucket_num + 1]; u++) { // TODO u umbenennen oder ohne ULLONG
        SETCHARBITPAIR(_g, u, 0);
    }

//    // Debug
//    cout << "g values for bucket " << bucket_num << " zeroed:" << endl;
//    for(ULLONG ii = _offset[bucket_num]; ii < _offset[bucket_num + 1]; ii++) {
//        cout << " " << GETCHARBITPAIR(_g, ii);
//    }
//    cout << endl;
//    // Debug end

    for(ULLONG j = next_queue_index - 1; j < next_queue_index; j--) {
        sum = 0;
        for(int k = 2; k >= 0; k--) {
            gValue = ARR(acyclicity_test_array, bucket_size, 3, queue[j], k);
            if(!GETBIT(visited, gValue)) {
                u = gValue;
                SETBIT(visited, gValue, 1);
                c = k;
            } else {
                sum += GETCHARBITPAIR(_g, _offset[bucket_num] + gValue);
            }
        }
        //SETCHARBITPAIR(_g, _offset[bucket_num] + u, (c - sum) % 3);
        SETCHARBITPAIR(_g, _offset[bucket_num] + u, (3 + c - (sum % 3)) % 3);
    }

//    // Debug
//    cout << "g values for bucket " << bucket_num << ":" << endl;
//    for(ULLONG ii = _offset[bucket_num]; ii < _offset[bucket_num + 1]; ii++) {
//        cout << " " << GETCHARBITPAIR(_g, ii);
//    }
//    cout << endl;
//    // Debug end

    delete[] edgesOf;
    delete[] cEdgesOf;
    delete[] queue;
    delete[] removed;
    delete[] visited;
    return false;
}

void PerfectHashFunction::_peelOf(ULLONG edge_index, ULLONG vertex_index, ULLONG *acyclicity_test_array,
                                  ULLONG bucket_size, ULLONG *queue, ULLONG &next_queue_index, ULLONG *edgesOf,
                                  ULLONG *cEdgesOf, ULLONG max_length, ULLONG mi, unsigned char *removed) {
    // TODO check this method
    ULLONG gValue, next_edge;

    //mark the edge as removed
    SETBIT(removed, edge_index, 1);
    queue[next_queue_index] = edge_index;
    next_queue_index++;

    //remove the edge from all incident nodes
    for(int k = 0; k < 3; k++) {
        gValue = ARR(acyclicity_test_array, bucket_size, 3, edge_index, k);
        if(gValue != vertex_index) {
            for(ULLONG l = 0; l < cEdgesOf[gValue]; l++) {
                if(ARR(edgesOf, mi, max_length, gValue, l) == edge_index) {
                    ARR(edgesOf, mi, max_length, gValue, l) = ARR(edgesOf, mi, max_length, gValue, cEdgesOf[gValue] - 1);
                    cEdgesOf[gValue]--;
                    break;
                }
            }

            //now do recursion on these nodes, if needed
            if(cEdgesOf[gValue] == 1) { //TODO: maybe we need to do this in a new loop
                next_edge = ARR(edgesOf, mi, max_length, gValue, 0);
                if(!GETBIT(removed, next_edge)) {
//                    // Debug
//                    cout << "remove edge " << next_edge << endl;
//                    // Debug end
                    _peelOf(next_edge, gValue, acyclicity_test_array, bucket_size, queue, next_queue_index, edgesOf,
                           cEdgesOf, max_length, mi, removed);
                }
            }
        }
    }
}

ULLONG PerfectHashFunction::evaluate(ULLONG x) {
    // TODO check this method
    ULLONG i, h0value, h1value, g0value, g1value, g2value, mi;
    ULLONG *h0coeffs = nullptr, *h1coeffs = nullptr;
    int sum;

    //find the bucket
    i = _evalUhf(x, _h_split_coeffs, _h_split_mod_mask, _m);
    h0coeffs = _h_coeffs + ((i * (_l + 1))<<1); // _h_coeffs + 2 * i * (_l + 1)
    h1coeffs = h0coeffs + _l + 1;
    mi = _offset[i+1] - _offset[i];

    //compute the values of hij(x)
    h0value = _evalUhf(x, h0coeffs, _h_mod_mask, _tab_rows);
    h1value = _evalUhf(x, h1coeffs, _h_mod_mask, _tab_rows);

    //compute the values of fij(x)
    //consider that we're adding the offset here as we need this later multiple times
    g0value = (((ARR(_random_table, _tab_rows, 6, h0value, 0) * _random_factor[i])
           ^ (ARR(_random_table, _tab_rows, 6, h1value, 1))) % mi) + _offset[i];
    g1value = (((ARR(_random_table, _tab_rows, 6, h0value, 2) * _random_factor[i])
           ^ (ARR(_random_table, _tab_rows, 6, h1value, 3))) % (mi - 1)) + _offset[i];
    g2value = (((ARR(_random_table, _tab_rows, 6, h0value, 4) * _random_factor[i])
           ^ (ARR(_random_table, _tab_rows, 6, h1value, 5))) % (mi - 2)) + _offset[i];

    //compute the values of the gij(x)
    if(g1value >= g0value) {
        g1value++;
    }
//    if(g2value >= g0value) {
//        if(g2value >= g1value) {
//            g2value += 2;
//        } else {
//            g2value++;
//        }
//    } else if(g2value >= g1value) {
//        g1value++;
//    }
    // first version not equivalent / not correct
    if(g2value >= g0value) {
        g2value++;
        if(g2value >= g1value) {
            g2value++;
        }
    } else if(g2value >= g1value) {
        g2value++;
        if(g2value >= g0value) {
            g2value++;
        }
    }

    //now compute the real hash value
    sum = (GETCHARBITPAIR(_g, g0value) + GETCHARBITPAIR(_g, g1value) + GETCHARBITPAIR(_g, g2value)) % 3;

    switch(sum) {
        case 0:
            return g0value;
        case 1:
            return g1value;
        default: //value 2
            return g2value;
    }
}

void PerfectHashFunction::_clear() {
    // TODO check this method
    delete[] _offset;
    delete[] _h_split_coeffs;
    delete[] _h_coeffs;
    delete[] _random_table;
    delete[] _random_factor;
    delete[] _g;
}

ULLONG PerfectHashFunction::getRange() {
    return _offset[_m];
}

ULLONG PerfectHashFunction::getSizeInBytes() {
    ULLONG size;
    size = 3 * sizeof(unsigned short) + 4 * sizeof(ULLONG);
    size += (_m + 1) * sizeof(ULLONG);
    size += (_l + 1) * sizeof(ULLONG);
    size += 2 * _m * (_l + 1) * sizeof(ULLONG);
    size += (ULLONG) ceil((long double)((6 * _tab_rows + _m) * _tab_width) / 8.0l);
    size += ((_offset[_m] >> 2) + 1) * sizeof(unsigned char);
    return size;
}