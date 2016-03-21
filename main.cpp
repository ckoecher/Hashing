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
#include "InputData.h"
#include "PerfectHashFunction.h"

/**
 * Creates a new data file if it doesn't exist yet and fills it with the first n integers.
 *
 * @param[in] dataFileName The path to the data file to create and fill
 */
void createInputData(string dataFileName) {
    cout << "\nCreate Input Data" << endl;
    cout << "n = ";
    ULLONG n;
    cin >> n;
    InputData *data = new InputData(dataFileName, ios::trunc);
    for (ULLONG i = 0; i < n; i++) {
        data->setNextValue(i);
    }
    data->close();
    delete data;
}

/**
 * Creates a new data file if it doesn't exist yet and fills it with random data.
 *
 * @param[in] dataFileName The path to the data file to create and fill
 */
void createRandomInputData(string dataFileName) {
    cout << "\nCreate Random Input Data" << endl;
    ULLONG n, seed;
    cout << "n = ";
    cin >> n;
    cout << "seed = ";
    cin >> seed;
    InputData *data = new InputData(dataFileName, ios::trunc);
    mt19937 *rng = new mt19937(seed);
    uniform_int_distribution<ULLONG> *dist = new uniform_int_distribution<ULLONG>;
    for (ULLONG i = 0; i < n; i++) {
        data->setNextValue((*dist)(*rng));
    }
    data->close();
    delete data;
    delete rng;
    delete dist;
}

/**
 * Reads the configuration file and inserts the data into a Configuration struct. Consider that the configuration file
 * has to be a INI file.
 *
 * @param[in] fileName The path to the configuration file
 * @return             The configuration data
 */
Configuration readConfig(string fileName) {
    INIReader reader(fileName);

    if (reader.ParseError() < 0) {
        cout << "Cannot load configuration file " << fileName << "." << endl << reader.ParseError() << endl;
        throw 0;
    }

    struct Configuration config;
    config.k = (unsigned short) reader.GetInteger("Hashing", "k", 32);
    config.l = (unsigned short) reader.GetInteger("Hashing", "l", 2);
    config.m_coeff = reader.GetReal("Hashing", "m_coeff", 2);
    config.m_exp = reader.GetReal("Hashing", "m_exp", 2.0 / 3);
    config.additional_bits_uhf = (unsigned short) reader.GetInteger("Hashing", "additional_bits_uhf", 6);
    config.num_of_tries_split = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_split", 42);
    config.num_of_tries_goodpairs = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_goodpairs", 42);
    config.mi_coeff = reader.GetReal("Hashing", "mi_coeff", 1.25);
    config.tab_rows_coeff = reader.GetReal("Hashing", "tab_rows_coeff", 2);
    config.tab_rows_exp = reader.GetReal("Hashing", "tab_rows_exp", 0.75);
    config.additional_bits_tab = (unsigned short) reader.GetInteger("Hashing", "additional_bits_tab", 6);
    config.num_of_tries_random_tab = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_random_tab", 42);
    config.num_of_tries_random_si = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_random_si", 42);
    config.seed = (ULLONG) reader.GetInteger("Hashing", "seed", 123456); //TODO which default value?
    config.debug_mode = reader.GetBoolean("Hashing", "debug_mode", true);

    return config;
}

/**
 * Checks whether the constructed perfect hash function is injective.
 *
 * @param[in,out] config The configuration data
 * @param[in]     data   The used data set
 * @param[in]     phf    The created perfect hash function
 * @param[in,out] stats  The additional statistical data
 * @return               true, if the function is injective, otherwise false
 */
bool validatePerfectHashFunction(Configuration &config, InputData *data, PerfectHashFunction *phf, Statistics &stats) {

    ULLONG datalength = data->getLength();
    ULLONG range = phf->getRange();
    ULLONG key;
    ULLONG hashvalue;
    ULLONG i;
    unsigned char *testarray = new unsigned char[(range >> 3) + 1]();
    bool perfect = true;

    // Stats
    stats.eval_start = clock();
    data->resetEvalTime();
    // Stats end

    // Debug
    if (config.debug_mode) {
        cout << "### Validate Perfect Hash Function ###" << endl;
        cout << "Check " << datalength << " keys:" << endl;
    }
    short percentage = -1;
    // Debug end

    for (i = 0; i < datalength; i++) {
        // Debug
        if (config.debug_mode) {
            if ((i * 100) / datalength != percentage) {
                percentage = (i * 100) / datalength;
                cout << "\r " << percentage << "% checked." << flush;
            }
        }
        // Debug end
        key = data->getValue(i);
        hashvalue = phf->evaluate(key);
        if (GETBIT(testarray, hashvalue)) {
            perfect = false;
            break;
        }
        SETBIT(testarray, hashvalue, 1);
    }
    // Stats
    stats.eval_end = clock();
    stats.eval_time = stats.eval_end - stats.eval_start;
    stats.eval_io = data->getEvalTime();
    // Stats end
    if (perfect) {
        // Debug
        if (config.debug_mode) {
            cout << "\r 100% checked." << endl;
            cout << "### Validate Perfect Hash Function Successful ###" << endl;
        }
        // Debug end
        // Stats
        stats.avg_eval_time = (long double) stats.eval_time / (long double) datalength;
        stats.avg_eval_io = (long double) stats.eval_io / (long double) datalength;
        stats.eval_success = true;
        // Stats end
        return true;
    } else {
        // Debug
        if (config.debug_mode) {
            cout << "\nKey no. " << i << " (" << key << ") with  hash value " << hashvalue << " collided." << endl;
            cout << "### Validate Perfect Hash Function Unsuccessful ###" << endl;
        }
        // Debug end
        return false;
    }
}

/**
 * Saves the statistical data of construction and validation of the perfect hash function into the statistics file.
 *
 * @param[in]     statsFileName  The path to the statistics file
 * @param[in]     configFileName The path to the configuration file
 * @param[in,out] config         The configuration data
 * @param[in]     dataFileName   The path to the data file
 * @param[in]     stats          The statistical data
 */
void saveStatistics(string statsFileName, string configFileName, Configuration &config, string dataFileName,
                    Statistics stats) {

    // (create and) open statistics file, new results will be appended
    fstream file(statsFileName, ios::out | ios::app);
    if(file.fail()) {
        cerr << "Could not open file \"" << statsFileName <<"\". Error: " << errno;
        throw 0;
    }

    // if file is empty, add header row
    if (file.tellg() == 0) {
        file << "0"
        << ";configFileName;m_coeff;m_exp;seed;dataFileName"
        << ";clocks_per_sec;num_of_keys;range_of_phf"
        << ";num_of_buckets;max_bucket_size;min_bucket_size;avg_bucket_size"
        << ";compact_size_in_bits;compact_size_overhead_in_bits;compact_size_overhead_percentaged"
        << ";compact_size_in_bits_general;compact_size_in_bits_split_uhf;compact_size_in_bits_offsets"
        << ";compact_size_in_bits_good_uhf_pairs;compact_size_in_bits_random_width;compact_size_in_bits_random_table"
        << ";compact_size_in_bits_random_factor;compact_size_in_bits_g_array"
        << ";creation_time;creation_io"
        << ";setup_time;setup_io"
        << ";split_time;split_io;split_tries;split_success"
        << ";goodpairs_time;goodpairs_io;goodpairs_total_tries;goodpairs_success"
        << ";buckets_time;buckets_io;random_tab_tries;random_si_total_tries;buckets_success"
        << ";eval_time;eval_io;avg_eval_time;avg_eval_io;eval_success"
        << ";size_in_bits"
        << ";size_in_bits_general;size_in_bits_split_uhf;size_in_bits_offsets"
        << ";size_in_bits_good_uhf_pairs;size_in_bits_random_width;size_in_bits_random_table"
        << ";size_in_bits_random_factor;size_in_bits_g_array"
        << endl;
    }

    // save statistics
    file << "1";
    file << ";\"" << configFileName << "\"";
    file << ";" << config.m_coeff;
    file << ";" << config.m_exp;
    file << ";" << config.seed;
    file << ";\"" << dataFileName << "\"";

    file << ";" << stats.clocks_per_sec;
    file << ";" << stats.num_of_keys;
    file << ";" << stats.range_of_phf;

    file << ";" << stats.num_of_buckets;
    file << ";" << stats.max_bucket_size;
    file << ";" << stats.min_bucket_size;
    file << ";" << stats.avg_bucket_size;

    file << ";" << stats.compact_size_in_bits;
    file << ";" << stats.compact_size_overhead_in_bits;
    file << ";" << stats.compact_size_overhead_percentaged;
    file << ";" << stats.compact_size_in_bits_general;
    file << ";" << stats.compact_size_in_bits_split_uhf;
    file << ";" << stats.compact_size_in_bits_offsets;
    file << ";" << stats.compact_size_in_bits_good_uhf_pairs;
    file << ";" << stats.compact_size_in_bits_random_width;
    file << ";" << stats.compact_size_in_bits_random_table;
    file << ";" << stats.compact_size_in_bits_random_factor;
    file << ";" << stats.compact_size_in_bits_g_array;

    //file << ";" << stats.creation_start;
    //file << ";" << stats.creation_end;
    file << ";" << stats.creation_time;
    file << ";" << stats.creation_io;
    //file << ";" << stats.creation_success;

    //file << ";" << stats.setup_start;
    //file << ";" << stats.setup_end;
    file << ";" << stats.setup_time;
    file << ";" << stats.setup_io;
    //file << ";" << stats.setup_success;

    //file << ";" << stats.split_start;
    //file << ";" << stats.split_end;
    file << ";" << stats.split_time;
    file << ";" << stats.split_io;
    file << ";" << stats.split_tries;
    file << ";" << stats.split_success;


    //file << ";" << stats.goodpairs_start;
    //file << ";" << stats.goodpairs_end;
    file << ";" << stats.goodpairs_time;
    file << ";" << stats.goodpairs_io;
    file << ";" << stats.goodpairs_total_tries;
    file << ";" << stats.goodpairs_success;

    //file << ";" << stats.buckets_start;
    //file << ";" << stats.buckets_end;
    file << ";" << stats.buckets_time;
    file << ";" << stats.buckets_io;
    file << ";" << stats.random_tab_tries;
    file << ";" << stats.random_si_total_tries;
    file << ";" << stats.buckets_success;

    //file << ";" << stats.eval_start;
    //file << ";" << stats.eval_end;
    file << ";" << stats.eval_time;
    file << ";" << stats.eval_io;
    file << ";" << stats.avg_eval_time;
    file << ";" << stats.avg_eval_io;
    file << ";" << stats.eval_success;

    file << ";" << stats.size_in_bits;
    file << ";" << stats.size_in_bits_general;
    file << ";" << stats.size_in_bits_split_uhf;
    file << ";" << stats.size_in_bits_offsets;
    file << ";" << stats.size_in_bits_good_uhf_pairs;
    file << ";" << stats.size_in_bits_random_width;
    file << ";" << stats.size_in_bits_random_table;
    file << ";" << stats.size_in_bits_random_factor;
    file << ";" << stats.size_in_bits_g_array;

    file << endl;
    file.close();
}

/**
 * Checks whether the file with the given path exists on file system yet. We need this dirty hack here are C++11 hasn't
 * a simple check for this.
 *
 * @param[in] name The path to the file to check
 * @return         true, if the file exists yet, otherwise false
 */
inline bool fileExists(const string name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

/**
 * The main function.
 *
 * @param[in] argc The number of arguments
 * @param[in] argv The arguments
 * @return         The return state
 */
int main(int argc, char *argv[]) {
    // Parse command line arguments
    // argv[0]: path and name of program itself
    // argv[1]: path of configuration file
    // argv[2]: path of data file
    // argv[3]: path of statistics file

    string configFileName;
    string dataFileName;
    string statsFileName;

    if (argc < 4) {
        cerr << "Too few command line arguments!" << endl;
        cerr << "./Hashing <configuration file> <data file> <statistics file>" << endl;
        return 1;
    }

    configFileName = argv[1];
    dataFileName = argv[2];
    statsFileName = argv[3];

    if(!fileExists(dataFileName)) {
        createRandomInputData(dataFileName);
    }

    // Setup
    Configuration config = readConfig(configFileName);
    // TODO some other conditions?
    if (config.k * config.l != sizeof(ULLONG) * 8) {
        cerr << "The parameters k and l have to satisfy the condition k * l == " << sizeof(ULLONG) * 8 << "." << endl;
        return 0;
    }
    InputData *data = new InputData(dataFileName);
    Statistics stats;

    // Creation and Validation
    PerfectHashFunction *phf;
    try {
        phf = new PerfectHashFunction(config, data, stats);
        validatePerfectHashFunction(config, data, phf, stats);
        delete phf;
    } catch (int e) {
        cerr << "Could not create Perfect Hash Function." << endl;
    }

    // Save statistics
    saveStatistics(statsFileName, configFileName, config, dataFileName, stats);

    // Cleanup
    delete data;

    return 0;
}