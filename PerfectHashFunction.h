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

#ifndef HASHING_PERFECTHASHFUNCTION_H
#define HASHING_PERFECTHASHFUNCTION_H

#include "definitions.h"
#include "InputData.h"

/**
 * This class is used for construction and evaluation of perfect hash functions using the split-and-share paradigm.
 */
class PerfectHashFunction {
private:
    // general information
    /**
     * The bit length of a single key fragment.
     */
    unsigned short _k;

    /**
     * The number of key fragments in a single key. Consider, that the length of a key is 64 bit. So this is 64 / _k.
     */
    unsigned short _l;

    // split information
    /**
     * The number of buckets the key set will be split into.
     */
    ULLONG _m;

    /**
     * The offsets of the first elements of each bucket in result data array. The size of a single bucket can be easily
     * computed from this array:
     *
     * @example m_i = _offset[i + 1] - _offset[i]
     *
     * Consider that this array contains _m + 1 values. Especially _offset[0] is 0 and _offset[_m] is the length of
     * input data stream multiplied with Configuration.mi_coeff.
     */
    ULLONG *_offset = nullptr;

    /**
     * A bit mask for better efficiency on evaluation of split function to avoid call of modulo operation.
     *
     * @see PerfectHashFunction::_evalUhf()
     */
    ULLONG _h_split_mod_mask;

    /**
     * The coefficients of the split function.
     */
    ULLONG *_h_split_coeffs = nullptr;

    // good 1-universal hash function pairs h^i_j (where 0<=i<_m and 0<=j<2) for each bucket information
    /**
     * A bit mask for better efficiency on evaluation of the hash function pairs h^i_j to avoid call of modulo operation.
     *
     * @see PerfectHashFunction::_evalUhf()
     */
    ULLONG _h_mod_mask;

    /**
     * The modulus of the hash function pairs h^i_j, i.e. these functions return values from range [_tab_rows].
     * Especially this is the number of rows of _random_table.
     *
     * @see PerfectHashFunction::_evalUhf()
     */
    ULLONG _tab_rows;

    /**
     * The coefficients of the hash function pairs h^i_j. The values are ordered in the following way:
     * h^0_0, h^0_1, h^1_0, h^1_1, h^2_0, h^2_1, ..., where "h^i_j" means an array of the coefficients of the hash
     * function h^i_j. So this array has length 2*_m*(_l+1) and contains values from range [_h_mod_mask + 1].
     */
    ULLONG *_h_coeffs = nullptr;

    // shared randomness information
    /**
     * The number of bits of each table entry.
     */
    unsigned short _tab_width;

    /**
     * The entries of the random tables used for computation of the hash functions f^i_j (where 0<=i<_m and
     * 0<=j<3). This is thought as a 2-dimensional array whose columns represent the single tables (in ordering
     * T^0_0, T^0_1, T^1_0, T^1_1, T^2_0, T^2_1). So it's an array of 6 * _tab_rows values of bit length _tab_width.
     */
    ULLONG *_random_table = nullptr;

    /**
     * The additional random factors for computation of the hash functions f^i_j. Since each bucket needs just one of
     * those values this array has length _m and contains values of bit length _tab_width.
     */
    ULLONG *_random_factor = nullptr;

    // selection arrays G_i where G_i: [m_i] -> [3], i in [_m]
    /**
     * The G_i arrays needed for computation of values of our perfect hash function. These are _offset[_m] many values
     * from range [3]. So we save 4 values into one char value. Those are ordered as follows:
     * G_0[0], G_0[1], G_0[2], ..., G_0[m_0-1], G_1[0], G_1[1], ..., G_(_m-1)[m_(_m-1)-1],
     * i.e. G_i[j] = _g[_offset[i]+j]. On evaluation we sum three values from this table to guarantee the injectivity of
     * the perfect hash function.
     */
    unsigned char *_g = nullptr;

    // debug information
    /**
     * This property indicates whether to output debug notes (true) or not (false).
     */
    bool _debug_mode;

    /**
     * Inserts the data from configuration file and computes the bit masks. This is step 1 in our construction plan.
     *
     * @param[in,out] config      The configuration data
     * @param[in]     data_length The length of the input file
     */
    void _configure(Configuration &config, ULLONG data_length);

    /**
     * Creates a (nearly) 1-universal hash function. This is step 2 in our construction plan.
     *
     * @param[out] coeffs The array of coefficients of the hash function
     * @param[in]  rng    The random number generator
     * @param[in]  dist   The distribution of the random numbers
     */
    void _createUhf(ULLONG *coeffs, mt19937 *rng, uniform_int_distribution<ULLONG> *dist);

    /**
     * Evaluates the 1-universal hash function. This is some kind of multiply-shift-hashing: At first we compute the
     * scalar product of the key and the coefficients, then we compute the residue of division by a power of 2 (we use
     * conjunction with the bit mask in third argument as this is more efficient). Afterwards we do a integer division
     * by another power of 2 (this is the known "multiply-shift-hashing"). Since this returns values from a range [2^k]
     * for some integer k (which probably isn't the range we need) we need to do another modulo operation with the
     * fourth argument as modulus.
     *
     * @param[in] key     The key from input data
     * @param[in] coeff   The coefficients of the hash function
     * @param[in] modMask The bit-mask of the inner modulo operation
     * @param[in] modulus The modulus of the outer modulo operation to fit into the correct range
     */
    ULLONG _evalUhf(ULLONG key, ULLONG *coeff, ULLONG modMask, ULLONG modulus);

    /**
     * Splits the input data with the help of the 1-universal hash function into buckets and checks whether these
     * buckets are small enough (that means smaller than the sqare root of input data length). These are the steps 4-7
     * in our construction plan. Additionally we compute _tab_width here.
     *
     * @param[in,out] config          The configuration data
     * @param[in]     data            The input data file
     * @param[in,out] bucket_data     The bucketed input data
     * @param[in,out] bucket_offsets  The offsets of the buckets in bucketed data
     * @param[out]    max_bucket_size The maximal size of the buckets
     * @param[out]    max_mi          The maximal range of the values of the second-level hash functions
     * @param[in,out] stats           The statistical data of the construction and evaluation process
     */
    bool _split(Configuration &config, InputData *data, InputData *bucket_data,
                ULLONG *&bucket_offsets, ULLONG *max_bucket_size, ULLONG *max_mi,
                Statistics &stats);

    /**
     * Creates a good pair of hash functions h^i_j (where 0<=i<_m and 0<=j<2) for each bucket. This is step 8 in our
     * construction plan. Consider that this method throws an exception if it's all tries to find good pairs failed.
     *
     * @param[in,out] config          The configuration data
     * @param[in]     bucket_data     The bucketed input data
     * @param[in]     bucket_offsets  The offsets of the buckets in bucketed data
     * @param[in]     max_bucket_size The maximal size of the buckets
     * @param[in]     rng             The random number generator
     * @param[in]     dist            The distribution of the random numbers
     * @param[in,out] stats           The statistical data of the construction and evaluation process
     */
    void _createGoodPairs(Configuration &config, InputData *bucket_data, ULLONG *bucket_offsets, ULLONG max_bucket_size,
                          mt19937 *rng, uniform_int_distribution<ULLONG> *dist, Statistics &stats);

    /**
     * Creates the random tables T^i_j, i.e. the array _random_table. This is step 9 in our construction plan.
     *
     * @param[in] rng  The random number generator
     * @param[in] dist The distribution of the random numbers
     */
    void _createRandomTables(mt19937 *rng, uniform_int_distribution<ULLONG> *dist);

    /**
     * Creates the random factor of a single bucket, i.e. _random_factor[bucket_num]. This is the first part of step 10
     * in our construction plan.
     *
     * @param[in] bucket_num The index of the bucket
     * @param[in] rng        The random number generator
     * @param[in] dist       The distribution of the random numbers
     */
    void _createRandomFactor(ULLONG bucket_num, mt19937 *rng, uniform_int_distribution<ULLONG> *dist);

    /**
     * Computes the hash functions f^i_j and g^i_j (where 0<=i<_m and 0<=j<3) from good hash function pair
     * (h^i_0, h^i_1), _random_table and _random_factor[bucket_num]. Additionally this method evaluations these
     * functions for all keys in the bucket. This is the other half of step 10 and the steps 11-12 in our construction
     * plan.
     *
     * @param[in]     bucket_num            The index of the bucket (bucket_num = i)
     * @param[in,out] acyclicity_test_array The image of the functions g^i_j of this bucket
     * @param[in]     bucket_data           The bucketed input data
     * @param[in]     bucket_offset         The offset of the first element in this bucket
     * @param[in]     bucket_size           The size of this bucket
     */
    void _computeGij(ULLONG bucket_num, ULLONG *acyclicity_test_array, InputData *bucket_data, ULLONG bucket_offset,
                     ULLONG bucket_size);

    /**
     * Computes the adjacency graph from the hyper-graph given in acyclicity_test_array and checks whether this graph is
     * cyclic. Additionally this method computes the array G_i. Afterwards this method resets the acyclicity_test_array
     * to zero, again. This is step 13 in out construction plan.
     *
     * @param[in]     bucket_num            The index of the bucket (bucket_num = i)
     * @param[in,out] acyclicity_test_array The image of the functions g^i_j of this bucket
     * @param[in]     bucket_size           The size of this bucket
     * @return                              true, if the adjacency graph is cyclic, otherwise false
     */
    bool _isCyclic(ULLONG bucket_num, ULLONG *acyclicity_test_array, ULLONG bucket_size);

    /**
     * This is a helper method for acyclicity check of the adjacency graph. This will remove an edge and an incident
     * node with degree 1 from this graph and call this method recursively on the other two incident nodes, if possible.
     *
     * @param[in]     edge_index            The edge to peel of here
     * @param[in]     vertex_index          The vertex with degree 1 to handle here
     * @param[in]     acyclicity_test_array The image of the functions g^i_j of this bucket
     * @param[in]     bucket_size           The size of this bucket
     * @param[in,out] queue                 The queue of the edges we've peeled of, yet
     * @param[in,out] next_queue_index      The next free index in the queue
     * @param[in,out] edgesOf               The map of nodes to their incident edges
     * @param[in,out] cEdgesOf              The map of nodes to the number of their incident edges
     * @param[in]     max_length            The maximal number of edges that a node is allowed to be incident to (very
     *                                      unlikely to be exceeded)
     * @param[in]     mi                    The number of nodes in this graph
     * @param[in,out] removed               An bitmap that indicates which edges are peeled of, yet
     */
    void _peelOf(ULLONG edge_index, ULLONG vertex_index, ULLONG *acyclicity_test_array, ULLONG bucket_size,
                 ULLONG *queue, ULLONG &next_queue_index, ULLONG *edgesOf, ULLONG *cEdgesOf, ULLONG max_length,
                 ULLONG mi, unsigned char *removed);

    /**
     * Computes the size of the actual memory we need for this perfect hash function. Additionally this method computes
     * the memory size of this perfect hash function if we would save values in a more space efficient way (e.g. in _g
     * array we use 2 bits for each element, even though we just need 1.5 bits for each element).
     *
     * @param[in,out] stats The statistical data of the construction and evaluation process
     */ //TODO: is this right???
    void _computeSizes(Statistics &stats);

    /**
     * Frees the memory which is used for this perfect hash function.
     */
    void _clear(); // delete data

public:
    /**
     * Constructor of this class. This method computes a new perfect hash function from input data.
     *
     * @param[in,out] config The configuration data
     * @param[in]     data   The input data
     * @param[in,out] stats  The statistical data of the construction and evaluation process
     */
    PerfectHashFunction(Configuration &config, InputData *data, Statistics &stats);

    /**
     * Evaluates the perfect hash function on the given key.
     *
     * @param[in] x The key
     * @return      The value
     */
    ULLONG evaluate(ULLONG x);

    /**
     * Returns the range of the image of this perfect hash function, i.e. this is
     * Configuration.mi_coeff * InputData::getLength() resp. _offset[_m].
     *
     * @return The range
     */
    ULLONG getRange();

    /**
     * Returns the size (in bytes) of the memory needed for this perfect hash function.
     *
     * @return The size
     */
    ULLONG getSizeInBytes();

    /**
     * Returns the size (in bytes) of the memory theoretically needed for this perfect hash function, if we would save
     * the data in a more space efficient way.
     *
     * @return The size
     */
    ULLONG getCompactSizeInBytes();

    /**
     * Destructor of this class. Just frees the memory of the arrays saved here.
     */
    virtual ~PerfectHashFunction() { _clear(); }
};

#endif //HASHING_PERFECTHASHFUNCTION_H