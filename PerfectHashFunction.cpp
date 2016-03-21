/*
 * Hashing - A simple implementation of Split-And-Share-Hashing.
 * Copyright (C) 2016  Philipp Schlag, Chris Köcher
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "definitions.h"
#include "PerfectHashFunction.h"

PerfectHashFunction::PerfectHashFunction(Configuration &config, InputData *data, Statistics &stats) {

    // Stats
    stats.creation_start = clock();
    stats.num_of_keys = data->getLength();
    stats.setup_start = clock();
    data->resetEvalTime();
    // Stats end

    // temporary variables and objects to compute the perfect hash function
    InputData *bucket_data = new InputData();   // temporary file to store the split data sorted in their buckets
    unsigned short num_of_tries_split = 0;      // current number of tries to split the data into small buckets
    bool split_successful;                      // indicates whether the split into buckets was successful (small buckets) or not
    ULLONG max_bucket_size;                     // maximal bucket size
    ULLONG max_mi;                              // maximal range of all (perfect) bucket hash functions
    ULLONG *bucket_offsets = nullptr;           // offsets of buckets to find correct block in the temporary file
    ULLONG *acyclicity_test_array = nullptr;    // array to represent a 3-graph and to test it for acyclicity
    bool badTables;                             // indicates whether the current table of random values has already been proven to be bad or not
    bool badFactor;                             // indicates whether the current random factor has already been proven to be bad or not
    unsigned short num_of_tries_tab = 0;        // current number of tries to find a good table of random values
    unsigned short num_of_tries_si;             // current number of tries to find a good random factor

    // random number generator (Mersenne twister) and distributions for random values within a specific range
    mt19937 *rng = nullptr;                                             // Mersenne twister (RNG)
    uniform_int_distribution<ULLONG> *dist_h_split_coeffs = nullptr;    // for split hash function
    uniform_int_distribution<ULLONG> *dist_h_coeffs = nullptr;          // for bucket hash function
    uniform_int_distribution<ULLONG> *dist_tables = nullptr;            // for table with random values and random factors

    // Debug
    if (config.debug_mode) {
        cout << "### Perfect Hash Function Creation ###" << endl;
        cout << data->getLength() << " key(s)." << endl;
        cout << "# Setup..." << endl;
    }
    // Debug end

    // read configuration and compute some parameters
    _configure(config, data->getLength());

    // create random number generator (Mersenne twister) with seed from configuration
    rng = new mt19937(config.seed);

    // Debug
    if (_debug_mode) {
        cout << "Setup successful." << endl;
        cout << "# Split data into buckets:" << endl;
    }
    // Debug end

    // Stats
    stats.setup_end = clock();
    stats.setup_time = stats.setup_end - stats.setup_start;
    stats.setup_io = data->getEvalTime();
    stats.setup_success = true;
    stats.split_start = clock();
    data->resetEvalTime();
    // Stats end

    // create distribution for random coefficients of the split hash function
    dist_h_split_coeffs = new uniform_int_distribution<ULLONG>(0, _h_split_mod_mask);

    // try to create a split hash function that creates small buckets
    _h_split_coeffs = new ULLONG[_l + 1];
    do {
        num_of_tries_split++;
        // guess an new split hash function
        _createUhf(_h_split_coeffs, rng, dist_h_split_coeffs);
        // try to split the set of keys into small buckets (in temporary file)
        split_successful = _split(config, data, bucket_data, bucket_offsets, &max_bucket_size, &max_mi, stats);
    } while (!split_successful && num_of_tries_split < config.num_of_tries_split);
    // TODO assert(ceil(log(max_mi)) <= 64)

    // abort if all splits were unsuccessful
    if (!split_successful) {
        // Debug
        if (_debug_mode) {
            cout << "\nCould not split data within " << num_of_tries_split << " tries." << endl;
        }
        // Debug end
        delete bucket_data;
        delete[] bucket_offsets;
        delete rng;
        delete dist_h_split_coeffs;
        _clear();
        throw 0;
    }

    // Debug
    if (_debug_mode) {
        cout << "\nCould split data within " << num_of_tries_split << " tries." << endl;
        cout << "Created " << _m << " buckets with <= " << max_bucket_size << " elements." << endl;
/*        cout << "Bucket sizes:" << endl;
        for (ULLONG i = 0; i < _m; i++) {
            cout << " " << bucket_offsets[i + 1] - bucket_offsets[i];
        }
        cout << endl;*/
    }
    // Debug end

    // Stats
    stats.split_end = clock();
    stats.split_time = stats.split_end - stats.split_start;
    stats.split_io = data->getEvalTime();
    stats.split_success = true;
    stats.split_tries = num_of_tries_split;
    stats.num_of_buckets = _m;
    stats.max_bucket_size = max_bucket_size;
    // stats.min_bucket_size computed in _split(...)
    stats.avg_bucket_size = (long double) stats.num_of_keys / (long double) stats.num_of_buckets;
    stats.goodpairs_start = clock();
    bucket_data->resetEvalTime();
    // Stats end

    // Debug
    if (_debug_mode) {
        cout << "# Create good pairs of universal hash functions:" << endl;
    }
    // Debug end

    // compute parameters for good pairs of 1-universal hash functions
    _tab_rows = (ULLONG) ceil(config.tab_rows_coeff * pow(max_bucket_size * max_bucket_size, config.tab_rows_exp));
    _h_mod_mask = (ULLONG) pow(2.0l, _k + ceil(log2(_tab_rows)) + config.additional_bits_uhf) - 1;

    // create distribution for random coefficients of the good hash function pairs
    dist_h_coeffs = new uniform_int_distribution<ULLONG>(0, _h_mod_mask);

    // try to create good hash function pairs for all buckets
    try {
        _createGoodPairs(config, bucket_data, bucket_offsets, max_bucket_size, rng, dist_h_coeffs, stats);
    } catch (int e) {
        // abort because creation was not possible
        delete bucket_data;
        delete[] bucket_offsets;
        delete rng;
        delete dist_h_split_coeffs;
        delete dist_h_coeffs;
        _clear();
        throw 0;
    }

    // Debug
    if (_debug_mode) {
        cout << "Created good pairs of universal hash functions for all " << _m << " buckets." << endl;
    }
    // Debug end

    // Stats
    stats.goodpairs_end = clock();
    stats.goodpairs_time = stats.goodpairs_end - stats.goodpairs_start;
    stats.goodpairs_io = bucket_data->getEvalTime();
    stats.goodpairs_success = true;
    stats.buckets_start = clock();
    bucket_data->resetEvalTime();
    // Stats end

    // Debug
    if (_debug_mode) {
        cout << "# Create perfect hash functions for all buckets:" << endl;
    }
    short percentage = -1;
    // Debug end

    // create distribution for shared random values (table and factors)
    _tab_width = (unsigned short) ceil(log2(max_bucket_size)) + config.additional_bits_tab;
    // TODO assert(_tab_width <= 64);
    dist_tables = new uniform_int_distribution<ULLONG>(0, (ULLONG) pow(2.0l, _tab_width) - 1);

    // prepare data structures to compute perfect bucket hash functions
    acyclicity_test_array = new ULLONG[max_mi * 3]();
    _random_table = new ULLONG[6 * _tab_rows];
    _random_factor = new ULLONG[_m];
    _g = new unsigned char[(_offset[_m] >> 2) + 1]();

    // try to create perfect hash functions for each bucket
    do {

        num_of_tries_tab++;
        // Stats
        stats.random_tab_tries++;
        // Stats end

        // create new table with random values
        _createRandomTables(rng, dist_tables);
        badTables = false;

        // try to create perfect hash functions for each bucket
        for (ULLONG i = 0; i < _m; i++) {
            // Debug
            if (_debug_mode) {
                if ((i * 100) / _m != percentage) {
                    percentage = (i * 100) / _m;
                    cout << "\r " << percentage << "% created." << flush;
                }
            }
            // Debug end
            num_of_tries_si = 0;

            // try to create a perfect hash function for bucket i
            do {
                num_of_tries_si++;
                // Stats
                stats.random_si_total_tries++;
                // Stats end
                // create a new random factor
                _createRandomFactor(i, rng, dist_tables);
                // create a 3-graph by computing 3 hash values for each key in the bucket
                _computeGij(i, acyclicity_test_array, bucket_data, bucket_offsets[i],
                            bucket_offsets[i + 1] - bucket_offsets[i]);
                // check if the created 3-graph is acyclic and transform the acyclic 3-graph into a perfect hash function for bucket i
                badFactor = _isCyclic(i, acyclicity_test_array, bucket_offsets[i + 1] - bucket_offsets[i]);
            } while (badFactor && num_of_tries_si < config.num_of_tries_random_si);

            // abort if all random factors were unsuccessful for bucket i, i.e. try a new table with random values
            if (badFactor) {
                badTables = true;
                break;
            }
        }
    } while (badTables && num_of_tries_tab < config.num_of_tries_random_tab);

    // Stats
    stats.buckets_end = clock();
    stats.buckets_time = stats.buckets_end - stats.buckets_start;
    stats.buckets_io = bucket_data->getEvalTime();
    stats.buckets_success = !badTables;
    // Stats end

    // cleanup of temporary data that is not needed anymore
    delete bucket_data;
    delete[] bucket_offsets;
    delete[] acyclicity_test_array;
    delete rng;
    delete dist_h_split_coeffs;
    delete dist_h_coeffs;
    delete dist_tables;

    // abort if all tables with random tables and all random factors could not create perfect hash functions for all buckets
    if (badTables) {
        // Debug
        if (_debug_mode) {
            cout << "\nCould not create perfect hash functions within " << num_of_tries_tab <<
            " tries (new random tables)." << endl;
        }
        // Debug end
        _clear();
        //TODO
        throw 0;
    }

    // Stats
    stats.creation_end = clock();
    stats.creation_time = stats.creation_end - stats.creation_start;
    stats.creation_io = stats.setup_io + stats.split_io + stats.goodpairs_io + stats.buckets_io;
    stats.creation_success = true;
    stats.range_of_phf = getRange();
    _computeSizes(stats);
    // Stats end

    // Debug
    if (_debug_mode) {
        cout << "\r 100% created." << endl;
        cout << "Created perfect hash functions for all " << _m << " buckets within " << num_of_tries_tab <<
        " tries (new random tables)." << endl;
        cout << "### Perfect Hash Function Creation Successful ###" << endl;
    }
    // Debug end
}

void PerfectHashFunction::_configure(Configuration &config, ULLONG data_length) {
    _k = config.k;
    _l = config.l;
    _m = min((ULLONG) ceil(config.m_coeff * pow(data_length, config.m_exp)),
             (ULLONG) ceil((long double) data_length / 20.0));
    _h_split_mod_mask = (ULLONG) pow(2.0l, _k + ceil(log2(_m)) + config.additional_bits_uhf) - 1;
    _debug_mode = config.debug_mode;
}

void PerfectHashFunction::_createUhf(ULLONG *coeffs, mt19937 *rng, uniform_int_distribution<ULLONG> *dist) {
    for (int i = 0; i < _l + 1; i++) {
        coeffs[i] = (*dist)(*rng);
    }
}

ULLONG PerfectHashFunction::_evalUhf(ULLONG key, ULLONG *coeff, ULLONG modMask, ULLONG modulus) {
    // & modMask == Mod (modMask+1) = Mod 2^c
    // modMask < 2^64 => x*y Mod 2^l == (x*y Mod 2^64) Mod 2^l
    // x*y Mod 2^64 == unsigned long long multiplication ("without" overflow)
    ULLONG keyMask = (ULLONG) pow(2.0l, _k) - 1;
    ULLONG res = coeff[_l] & modMask;
    for (int i = _l - 1; i >= 0; i--) {
        res += (coeff[i] * (key & keyMask)) & modMask;
        res &= modMask;
        key >>= _k;
    }
    res >>= _k;
    res %= modulus;
    return res;
}

bool PerfectHashFunction::_split(Configuration &config, InputData *data, InputData *bucket_data,
                                 ULLONG *&bucket_offsets, ULLONG *max_bucket_size, ULLONG *max_mi, Statistics &stats) {
    // temporary variables and arrays to split the set of keys into small buckets
    ULLONG data_length = data->getLength();         // number of keys
    ULLONG *bucket_sizes = new ULLONG[_m]();        // bucket_sizes[i] == number of keys sorted into bucket i by h_split
    ULLONG hv;                                      // some hash value h_split(key)
    ULLONG mi_1;                                    // for computation of max_mi and _offset
    ULLONG bucketOverflowSize = max((ULLONG) floor(sqrt(data_length)), (ULLONG) 40);   // size of a bucket that is not small anymore
    ULLONG value;                                   // temporary value from data stream (key)

    // Debug
    short percentage = -1;
    // Debug end

    // counting
    for (ULLONG i = 0; i < data_length; i++) {
        // Debug
        if (_debug_mode) {
            if ((i * 100) / data_length != percentage) {
                percentage = (i * 100) / data_length;
                cout << "\r " << percentage << "% counted." << flush;
            }
        }
        // Debug end
        hv = _evalUhf(data->getValue(i), _h_split_coeffs, _h_split_mod_mask, _m);
        if (bucket_sizes[hv] >= bucketOverflowSize) {
            // bucket i overflow (too large)
            delete[] bucket_sizes;
            return false;
        } else {
            bucket_sizes[hv]++;
        }
    }

    // Debug
    if (_debug_mode) {
        cout << "\r 100% counted." << flush;
        cout << endl;
    }
    // Debug end

    // compute maximal bucket size and offsets (for next segments) for sorted data (bucket_data segments)
    bucket_offsets = new ULLONG[_m + 1]();
    bucket_offsets[0] = bucket_sizes[0];
    *max_bucket_size = bucket_sizes[0];
    // Stats
    stats.min_bucket_size = bucket_sizes[0];
    // Stats end
    for (int i = 1; i < _m; i++) {
        bucket_offsets[i] = bucket_offsets[i - 1] + bucket_sizes[i];
        if (*max_bucket_size < bucket_sizes[i]) {
            *max_bucket_size = bucket_sizes[i];
        }
        // Stats
        if (stats.min_bucket_size > bucket_sizes[i]) {
            stats.min_bucket_size = bucket_sizes[i];
        }
        // Stats end
    }
    bucket_offsets[_m] = bucket_offsets[_m - 1];

    // Debug
    percentage = -1;
    // Debug end

    // sort data into buckets (in temporary file)
    for (ULLONG i = data_length - 1; i < data_length; i--) {
        // Debug
        if (_debug_mode) {
            if (((data_length - 1 - i) * 100) / data_length != percentage) {
                percentage = ((data_length - 1 - i) * 100) / data_length;
                cout << "\r " << percentage << "% split." << flush;
            }
        }
        // Debug end
        value = data->getValue(i);
        hv = _evalUhf(value, _h_split_coeffs, _h_split_mod_mask, _m);
        bucket_data->setValue(value, bucket_offsets[hv] - 1);
        bucket_offsets[hv]--;
    }

    // Debug
    if (_debug_mode) {
        cout << "\r 100% split." << flush;
    }
    // Debug end

    // compute offsets for ranges of bucket hash functions and maximal range of a bucket hash function
    _offset = new ULLONG[_m + 1];
    _offset[0] = 0;
    *max_mi = 0;
    for (ULLONG i = 1; i < _m + 1; i++) {
        mi_1 = (ULLONG) ceil(config.mi_coeff * bucket_sizes[i - 1]);
        // to eliminate problem with %(mi-2) and
        // to eliminate problem with too small buckets and acyclicity
        if (mi_1 > 0 && mi_1 < 10) {
            mi_1 = 10;
        }
        if (*max_mi < mi_1) {
            *max_mi = mi_1;
        }
        _offset[i] = _offset[i - 1] + mi_1;
    }

    delete[] bucket_sizes;
    return true;
}

void PerfectHashFunction::_createGoodPairs(Configuration &config, InputData *bucket_data, ULLONG *bucket_offsets,
                                           ULLONG max_bucket_size, mt19937 *rng, uniform_int_distribution<ULLONG> *dist,
                                           Statistics &stats) {
    ULLONG *hashValues = new ULLONG[max_bucket_size << 1];              // array to store hashvalues of keys created by good pair of 1-universal hash function of current bucket
    unsigned char *hTables = new unsigned char[(_tab_rows >> 1) + 1](); // array to count number of keys with the same hash value (by the same hash function)
    ULLONG *h0coeffs = nullptr;                                         // coefficients of first hash function of current good pair of hash functions
    ULLONG *h1coeffs = nullptr;                                         // coefficients of second hash function of current good pair of hash functions
    bool goodPair;                                                      // indicates whether the current pair of hash functions is still considered as good
    ULLONG x;                                                           // current key
    ULLONG bucket_size;                                                 // size of current bucket
    unsigned short num_of_tries_goodpairs;                              // current number of tries to find a good pair (for each bucket)

    // Debug
    short percentage = -1;
    // Debug end

    // try to find a good pair of 1-universal hash functions for each bucket
    _h_coeffs = new ULLONG[(_m * (_l + 1)) << 1];
    for (ULLONG pairI = 0; pairI < _m; pairI++) {

        num_of_tries_goodpairs = 0;

        // find location to store the coefficients of the good pair of hash functions
        h0coeffs = _h_coeffs + ((pairI * (_l + 1)) << 1);
        h1coeffs = h0coeffs + _l + 1;

        // size of bucket pairI
        bucket_size = bucket_offsets[pairI + 1] - bucket_offsets[pairI];

        // Debug
        if (_debug_mode) {
            if ((pairI * 100) / _m != percentage) {
                percentage = (pairI * 100) / _m;
                cout << "\r " << percentage << "% created." << flush;
            }
        }
        // Debug end

        // try to find a good pair of 1-universal hash functions for bucket pairI
        do {
            goodPair = true;
            num_of_tries_goodpairs++;
            // Stats
            stats.goodpairs_total_tries++;
            // Stats end

            // guess an new pair of hash function
            _createUhf(h0coeffs, rng, dist);
            _createUhf(h1coeffs, rng, dist);

            // counting (for each hash function separately): how many keys have a specific hashvalue
            for (ULLONG j = 0; j < bucket_size; j++) {
                // compute hashvalues
                x = bucket_data->getValue(bucket_offsets[pairI] + j);
                ARR(hashValues, bucket_size, 2, j, 0) = _evalUhf(x, h0coeffs, _h_mod_mask, _tab_rows);
                ARR(hashValues, bucket_size, 2, j, 1) = _evalUhf(x, h1coeffs, _h_mod_mask, _tab_rows);
                // count hashvalues
                if (GETTABBITPAIR(hTables, 0, ARR(hashValues, bucket_size, 2, j, 0)) < 2) {
                    INCTABBITPAIR(hTables, 0, ARR(hashValues, bucket_size, 2, j, 0));
                }
                if (GETTABBITPAIR(hTables, 1, ARR(hashValues, bucket_size, 2, j, 1)) < 2) {
                    INCTABBITPAIR(hTables, 1, ARR(hashValues, bucket_size, 2, j, 1));
                }
            }
            // evaluation: does every key has at least one unique hashvalue by a specific hash function?
            for (ULLONG j = 0; j < bucket_size && goodPair; j++) {
                // TODO &&
                if ((GETTABBITPAIR(hTables, 0, ARR(hashValues, bucket_size, 2, j, 0)) == 2)
                    && (GETTABBITPAIR(hTables, 1, ARR(hashValues, bucket_size, 2, j, 1)) == 2)) {
                    // pair of hash functions is not good
                    goodPair = false;
                }
            }
            // clear hTables for new counting
            // TODO without condition?
            if (pairI < _m - 1 || !goodPair) {
                for (ULLONG j = 0; j < bucket_size; j++) {
                    ZEROTABBITPAIRS(hTables, ARR(hashValues, bucket_size, 2, j, 0));
                    ZEROTABBITPAIRS(hTables, ARR(hashValues, bucket_size, 2, j, 1));
                }
            }
        } while (!goodPair && num_of_tries_goodpairs < config.num_of_tries_goodpairs);

        // no good pair of hash function found for bucket pairI
        if (!goodPair) {
            // Debug
            if (_debug_mode) {
                cout << "\nCould not create good pair of universal hash functions for bucket " << pairI
                << " within " << num_of_tries_goodpairs << " tries." << endl;
            }
            // Debug end
            delete[] hTables;
            delete[] hashValues;
            // TODO
            throw 0;
        }
    }
    // Debug
    if (_debug_mode) {
        cout << "\r 100% created." << endl;
    }
    // Debug end
    delete[] hTables;
    delete[] hashValues;
}

void PerfectHashFunction::_createRandomTables(mt19937 *rng, uniform_int_distribution<ULLONG> *dist) {
    for (int i = 0; i < 6 * _tab_rows; i++) {
        _random_table[i] = (*dist)(*rng);
    }
}

void PerfectHashFunction::_createRandomFactor(ULLONG bucket_num, mt19937 *rng, uniform_int_distribution<ULLONG> *dist) {
    _random_factor[bucket_num] = (*dist)(*rng);
}

void PerfectHashFunction::_computeGij(ULLONG bucket_num, ULLONG *acyclicity_test_array,
                                      InputData *bucket_data, ULLONG bucket_offset, ULLONG bucket_size) {
    ULLONG h0value; // hash value of first hash function of current good pair of hash functions
    ULLONG h1value; // hash value of second hash function of current good pair of hash functions
    ULLONG fi0;     // first possible hash value for (perfect) bucket hash function
    ULLONG fi1;     // second possible hash value for (perfect) bucket hash function
    ULLONG fi2;     // third possible hash value for (perfect) bucket hash function
    ULLONG *h0coeffs = _h_coeffs + ((bucket_num * (_l + 1)) << 1);  // coefficients of first hash function of current good pair of hash functions
    ULLONG *h1coeffs = h0coeffs + _l + 1;                           // coefficients of second hash function of current good pair of hash functions
    ULLONG x;       // current key
    ULLONG mi = _offset[bucket_num + 1] - _offset[bucket_num];      // range of (perfect) bucket hash function

    // compute three hash values for each key and construct a 3-graph (stored in acyclicity_test_array)
    for (ULLONG k = 0; k < bucket_size; k++) {
        // compute hash values of good pair of 1-universal hash functions
        x = bucket_data->getValue(bucket_offset + k);
        h0value = _evalUhf(x, h0coeffs, _h_mod_mask, _tab_rows);
        h1value = _evalUhf(x, h1coeffs, _h_mod_mask, _tab_rows);

        // compute the values fij(x)
        fi0 = ((ARR(_random_table, _tab_rows, 6, h0value, 0) * _random_factor[bucket_num])
               ^ (ARR(_random_table, _tab_rows, 6, h1value, 1))) % mi;
        fi1 = ((ARR(_random_table, _tab_rows, 6, h0value, 2) * _random_factor[bucket_num])
               ^ (ARR(_random_table, _tab_rows, 6, h1value, 3))) % (mi - 1);
        fi2 = ((ARR(_random_table, _tab_rows, 6, h0value, 4) * _random_factor[bucket_num])
               ^ (ARR(_random_table, _tab_rows, 6, h1value, 5))) % (mi - 2);

        // compute the pairwise distinct values gij(x) as possible results of the (perfect) bucket hash function
        if (fi1 >= fi0) {
            fi1++;
        }
        if (fi2 >= fi0) {
            fi2++;
            if (fi2 >= fi1) {
                fi2++;
            }
        } else if (fi2 >= fi1) {
            fi2++;
            if (fi2 >= fi0) {
                fi2++;
            }
        }

        // construct the 3-graph by adding the edge (fi0, fi1, fi2)
        ARR(acyclicity_test_array, bucket_size, 3, k, 0) = fi0;
        ARR(acyclicity_test_array, bucket_size, 3, k, 1) = fi1;
        ARR(acyclicity_test_array, bucket_size, 3, k, 2) = fi2;
    }
}

bool PerfectHashFunction::_isCyclic(ULLONG bucket_num, ULLONG *acyclicity_test_array, ULLONG bucket_size) {
    ULLONG max_length = 2 * (ULLONG) ceil(log2(bucket_size)) + 1;   // maximal number of edges that a node is allowed to be incident to (very unlikely to be exceeded)
    ULLONG mi = _offset[bucket_num + 1] - _offset[bucket_num];      // number of nodes of the 3-graph = range of the bucket hash function
    ULLONG *edgesOf = new ULLONG[max_length * mi]();                // array to store the incident edges of each node
    ULLONG *cEdgesOf = new ULLONG[mi]();                            // array to count the number of incident edges of each node
    ULLONG gValue;                                                  // node of the 3-graph = (possible) hashvalue of the bucket hash function
    ULLONG edge_index;                                              // index of an edge of the 3-graph = index of a key of the bucket
    ULLONG *queue = nullptr;                                        // stores all edges (keys) that have been removed from the 3-graph in order to test for acyclicity
    ULLONG next_queue_index = 0;                                    // number of already removed edges
    unsigned char *removed = nullptr;                               // bitmap that stores information if an edge has been removed
    unsigned char *visited = nullptr;                               // bitmap that stores information if a nodes has been visited already
    ULLONG u;                                                       // stores a newly visited node in the construction phase of the _g array
    int c, sum;                                                     // variables to compute the correct value for the _g array for each key (edge)

    // construct a list of incident edges of each node
    for (ULLONG j = 0; j < bucket_size; j++) {
        for (int k = 0; k < 3; k++) {
            gValue = ARR(acyclicity_test_array, bucket_size, 3, j, k);
            // abort if node gValue has too many incident edges
            if (cEdgesOf[gValue] >= max_length) {
                delete[] edgesOf;
                delete[] cEdgesOf;
                return true;
            }
            ARR(edgesOf, mi, max_length, gValue, cEdgesOf[gValue]) = j;
            cEdgesOf[gValue]++;
        }
    }

    // check for acyclicity
    queue = new ULLONG[bucket_size];
    removed = new unsigned char[(bucket_size >> 3) + 1]();
    // check for each node (hashvalue)
    for (ULLONG j = 0; j < mi; j++) {
        // if the current node has only one remaining incident edge (only one remaining key has this hashvalue)
        if (cEdgesOf[j] == 1) {
            // remove the edge (key) (if possible) and check the nodes (hashvalues) (formerly) incident to the removed edge
            edge_index = ARR(edgesOf, mi, max_length, j, 0);
            if (!GETBIT(removed, edge_index)) {
                _peelOf(edge_index, j, acyclicity_test_array, bucket_size, queue, next_queue_index, edgesOf,
                        cEdgesOf, max_length, mi, removed);
            }
        }
    }

    // abort if there are still some edges (keys) that have not been removed (incident nodes are incident to at least one other edge)
    // i.e. abort if the 3-graph is not acyclic
    if (next_queue_index != bucket_size) {
        delete[] edgesOf;
        delete[] cEdgesOf;
        delete[] queue;
        delete[] removed;
        return true;
    }

    // now assign the necessary values to _g in order to transform the acyclic 3-graph into a perfect hash function
    visited = new unsigned char[(mi >> 3) + 1]();
    for (ULLONG v = _offset[bucket_num]; v < _offset[bucket_num + 1]; v++) {
        SETCHARBITPAIR(_g, v, 0);
    }
    for (ULLONG j = next_queue_index - 1; j < next_queue_index; j--) {
        sum = 0;
        for (int k = 2; k >= 0; k--) {
            gValue = ARR(acyclicity_test_array, bucket_size, 3, queue[j], k);
            if (!GETBIT(visited, gValue)) {
                u = gValue;
                SETBIT(visited, gValue, 1);
                c = k;
            } else {
                sum += GETCHARBITPAIR(_g, _offset[bucket_num] + gValue);
            }
        }
        SETCHARBITPAIR(_g, _offset[bucket_num] + u, (3 + c - (sum % 3)) % 3);
    }

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
    ULLONG gValue;      // node of the 3-graph = (possible) hashvalue of the bucket hash function
    ULLONG next_edge;   // index of an edge of the 3-graph = index of a key of the bucket

    // mark the edge as removed
    SETBIT(removed, edge_index, 1);
    queue[next_queue_index] = edge_index;
    next_queue_index++;

    // remove the edge from all incident nodes and do recursion in these nodes (if needed)
    for (int k = 0; k < 3; k++) {
        // node component number k (of 3) of the edge
        gValue = ARR(acyclicity_test_array, bucket_size, 3, edge_index, k);
        // check if this node component is not the node from the function call
        if (gValue != vertex_index) {
            // remove the edge from the list of incident edges
            for (ULLONG l = 0; l < cEdgesOf[gValue]; l++) {
                if (ARR(edgesOf, mi, max_length, gValue, l) == edge_index) {
                    ARR(edgesOf, mi, max_length, gValue, l) = ARR(edgesOf, mi, max_length, gValue,
                                                                  cEdgesOf[gValue] - 1);
                    cEdgesOf[gValue]--;
                    break;
                }
            }
            // now do recursion on this node if he has only one remaining incident edge
            if (cEdgesOf[gValue] == 1) {
                next_edge = ARR(edgesOf, mi, max_length, gValue, 0);
                if (!GETBIT(removed, next_edge)) {
                    _peelOf(next_edge, gValue, acyclicity_test_array, bucket_size, queue, next_queue_index, edgesOf,
                            cEdgesOf, max_length, mi, removed);
                }
            }
        }
    }
}

void PerfectHashFunction::_computeSizes(Statistics &stats) {
    stats.size_in_bits_general = 8* 2 * sizeof(unsigned short);
    stats.size_in_bits_split_uhf = 8 * (_l + 3) * sizeof(ULLONG);
    stats.size_in_bits_offsets = 8 * (_m + 1) * sizeof(ULLONG);
    stats.size_in_bits_good_uhf_pairs = 8 * (2 * _m * (_l + 1) + 2) * sizeof(ULLONG);
    stats.size_in_bits_random_width = 8 * sizeof(unsigned short);
    stats.size_in_bits_random_table = 8 * 6 * _tab_rows * sizeof(ULLONG);
    stats.size_in_bits_random_factor = 8 * _m * sizeof(ULLONG);
    stats.size_in_bits_g_array = 8 * ((_offset[_m] >> 2) + 1) * sizeof(unsigned char);
    stats.size_in_bits = stats.size_in_bits_general + stats.size_in_bits_offsets + stats.size_in_bits_split_uhf
                          + stats.size_in_bits_good_uhf_pairs + stats.size_in_bits_random_width
                          + stats.size_in_bits_random_table + stats.size_in_bits_random_factor
                          + stats.size_in_bits_g_array;

    stats.compact_size_in_bits_general = 8 * 2 * sizeof(unsigned short);
    stats.compact_size_in_bits_split_uhf =
            8 * 2 * sizeof(ULLONG) + (_l + 1) * (ULLONG)log2(_h_split_mod_mask + 1);
    stats.compact_size_in_bits_offsets = 8 * (_m + 1) * sizeof(ULLONG);
    stats.compact_size_in_bits_good_uhf_pairs =
            8 * 2 * sizeof(ULLONG) + (2 * _m * (_l + 1)) * (ULLONG)log2(_h_mod_mask + 1);
    stats.compact_size_in_bits_random_width = 8 * sizeof(unsigned short);
    stats.compact_size_in_bits_random_table = 6 * _tab_rows * _tab_width;
    stats.compact_size_in_bits_random_factor = _m * _tab_width;
    stats.compact_size_in_bits_g_array = 8 * ((_offset[_m] >> 2) + 1) * sizeof(unsigned char);
    stats.compact_size_in_bits = stats.compact_size_in_bits_general + stats.compact_size_in_bits_offsets
                                  + stats.compact_size_in_bits_split_uhf + stats.compact_size_in_bits_good_uhf_pairs
                                  + stats.compact_size_in_bits_random_width + stats.compact_size_in_bits_random_table
                                  + stats.compact_size_in_bits_random_factor + stats.compact_size_in_bits_g_array;
    stats.compact_size_overhead_in_bits = stats.compact_size_in_bits - (ULLONG)ceil(2.5 * stats.num_of_keys);
    stats.compact_size_overhead_percentaged = 100 * (long double)stats.compact_size_overhead_in_bits
                                              / (long double)(2.5 * stats.num_of_keys);
}

void PerfectHashFunction::_clear() {
    delete[] _offset;
    delete[] _h_split_coeffs;
    delete[] _h_coeffs;
    delete[] _random_table;
    delete[] _random_factor;
    delete[] _g;
}

ULLONG PerfectHashFunction::evaluate(ULLONG x) {
    ULLONG i;                   // bucket number for key x
    ULLONG mi;                  // offset of range of bucket i
    ULLONG *h0coeffs = nullptr; // coefficients of first hash function of good pair of hash functions of bucket i
    ULLONG *h1coeffs = nullptr; // coefficients of second hash function of good pair of hash functions of bucket i
    ULLONG h0value;             // hash value of first hash function of good pair of hash functions of bucket i
    ULLONG h1value;             // hash value of second hash function of good pair of hash functions of bucket i
    ULLONG g0value;             // first possible hash value of perfect bucket hash function
    ULLONG g1value;             // second possible hash value of perfect bucket hash function
    ULLONG g2value;             // third possible hash value of perfect bucket hash function
    int sum;                    // variable to determine the correct hash value of the perfect bucket hash function

    // find the bucket, its good pair of hash function and the offset of its range
    i = _evalUhf(x, _h_split_coeffs, _h_split_mod_mask, _m);
    h0coeffs = _h_coeffs + ((i * (_l + 1)) << 1); // _h_coeffs + 2 * i * (_l + 1)
    h1coeffs = h0coeffs + _l + 1;
    mi = _offset[i + 1] - _offset[i];

    // compute the values hij(x)
    h0value = _evalUhf(x, h0coeffs, _h_mod_mask, _tab_rows);
    h1value = _evalUhf(x, h1coeffs, _h_mod_mask, _tab_rows);

    // compute the values fij(x)
    // consider that we're adding the offset here as we need this later multiple times
    g0value = (((ARR(_random_table, _tab_rows, 6, h0value, 0) * _random_factor[i])
                ^ (ARR(_random_table, _tab_rows, 6, h1value, 1))) % mi) + _offset[i];
    g1value = (((ARR(_random_table, _tab_rows, 6, h0value, 2) * _random_factor[i])
                ^ (ARR(_random_table, _tab_rows, 6, h1value, 3))) % (mi - 1)) + _offset[i];
    g2value = (((ARR(_random_table, _tab_rows, 6, h0value, 4) * _random_factor[i])
                ^ (ARR(_random_table, _tab_rows, 6, h1value, 5))) % (mi - 2)) + _offset[i];

    // compute the values gij(x)
    if (g1value >= g0value) {
        g1value++;
    }
    if (g2value >= g0value) {
        g2value++;
        if (g2value >= g1value) {
            g2value++;
        }
    } else if (g2value >= g1value) {
        g2value++;
        if (g2value >= g0value) {
            g2value++;
        }
    }

    // compute the real hash value
    sum = (GETCHARBITPAIR(_g, g0value) + GETCHARBITPAIR(_g, g1value) + GETCHARBITPAIR(_g, g2value)) % 3;
    switch (sum) {
        case 0:
            return g0value;
        case 1:
            return g1value;
        default: // sum == 2
            return g2value;
    }
}

ULLONG PerfectHashFunction::getRange() {
    return _offset[_m];
}

ULLONG PerfectHashFunction::getSizeInBits() {
    ULLONG size;
    size = 8 * (3 * sizeof(unsigned short) + 4 * sizeof(ULLONG));
    size += 8 * (_m + 1) * sizeof(ULLONG);
    size += 8 * (_l + 1) * sizeof(ULLONG);
    size += 8 * 2 * _m * (_l + 1) * sizeof(ULLONG);
    size += 8 * (6 * _tab_rows + _m) * sizeof(ULLONG);
    size += 8 * ((_offset[_m] >> 2) + 1) * sizeof(unsigned char);
    return size;
}

ULLONG PerfectHashFunction::getCompactSizeInBits() {
    ULLONG size;
    size = 8 * (3 * sizeof(unsigned short) + 4 * sizeof(ULLONG));
    size += 8 * (_m + 1) * sizeof(ULLONG);
    size += (_l + 1) * (ULLONG)log2(_h_split_mod_mask + 1);
    size += (2 * _m * (_l + 1)) * (ULLONG)log2(_h_mod_mask + 1);
    size += (6 * _tab_rows + _m) * _tab_width;
    size += 8 * ((_offset[_m] >> 2) + 1) * sizeof(unsigned char);
    return size;
}