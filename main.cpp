#include "definitions.h"
#include "InputData.h"
#include "PerfectHashFunction.h"

void testTypeSizes() {
    cout << "\nBytes of unsigned char: " << sizeof(unsigned char) << endl;
    assert(sizeof(unsigned char) == 1);
    cout << "Bytes of ULLONG: " << sizeof(ULLONG) << endl;
    assert(sizeof(ULLONG) == 8);
    cout << "Bytes of double: " << sizeof(double) << ", long double: " << sizeof(long double) << endl;
    assert(sizeof(long double) == 16);
}

void testMersenneTwister() {
    unsigned long seed;
    cout << "\nseed = ";
    cin >> seed;
    mt19937 *rng;
    rng = new mt19937(seed);
    uniform_int_distribution<ULLONG> *dist_test;
    dist_test = new uniform_int_distribution<ULLONG>(0, 1);
    for (int i = 0; i < 10; i++) {
        cout << (*dist_test)(*rng) << endl;
    }
    delete rng;
    delete dist_test;
}

void testBitPairs() {
    cout << "\nTest GETTABBITPAIR, ZEROTABBITPAIRS and INCTABBITPAIR: carray[16]";
    int n = 16;
    clock_t cl;
    unsigned char *carray = new unsigned char[n]();
    unsigned char c;
    cout << "carray:";
    for (int i = 0; i < n; i++) {
        carray[i] = i;
        cout << carray[i];
    }
    cout << "\nGETTABBITPAIR content:";
    for (int i = 0; i < 2 * n; i++) {
        for (int j = 0; j < 2; j++) {
            cout << " " << GETTABBITPAIR(carray, j, i);
        }
        cout << " #";
    }
    cout << "\nZEROTABBITPAIRS content:";
    for (int i = 0; i < 2 * n; i++) {
        for (int j = 0; j < 2; j++) {
            ZEROTABBITPAIRS(carray, i);
        }
    }
    for (int i = 0; i < 2 * n; i++) {
        for (int j = 0; j < 2; j++) {
            cout << " " << GETTABBITPAIR(carray, j, i);
        }
        cout << " #";
    }
    cout << "\nINCTABBITPAIR content:";
    for (int i = 0; i < 2 * n; i++) {
        for (int j = 0; j < 2; j++) {
            if (i % 2 == 1) {
                INCTABBITPAIR(carray, j, i);
                INCTABBITPAIR(carray, j, i);
            }
            if (j == 1) {
                INCTABBITPAIR(carray, j, i);
            }
        }
    }
    for (int i = 0; i < 2 * n; i++) {
        for (int j = 0; j < 2; j++) {
            cout << " " << GETTABBITPAIR(carray, j, i);
        }
        cout << " #";
    }
    cout << "\nn = ";
    cin >> n;
    delete[] carray;
    carray = new unsigned char[n];
    for (int round = 0; round < 5; round++) {
        cout << "\nGETTABBITPAIR time to read: ";
        cl = clock();
        for (int i = 0; i < 2 * n; i++) {
            for (int j = 0; j < 2; j++) {
                c = GETTABBITPAIR(carray, j, i);
            }
        }
        cout << clock() - cl << " ticks";
    }
    delete[] carray;
}

void testSetGetBit() {
    cout << "\nTest SETBIT and GETBIT" << endl;
    unsigned char *array = new unsigned char[2]();
    for (int i = 0; i < 2 * sizeof(unsigned char) * 8; i++) {
        cout << GETBIT(array, i) << " ";
    }
    cout << endl;
    for (int j = 0; j < 2 * sizeof(unsigned char) * 8; j++) {
        SETBIT(array, j, 1);
        for (int i = 0; i < 2 * sizeof(unsigned char) * 8; i++) {
            cout << GETBIT(array, i) << " ";
        }
        cout << endl;
    }
    for (int j = 0; j < 2 * sizeof(unsigned char) * 8; j++) {
        SETBIT(array, j, 0);
        for (int i = 0; i < 2 * sizeof(unsigned char) * 8; i++) {
            cout << GETBIT(array, i) << " ";
        }
        cout << endl;
    }
}

void testSetGetCharBitPairs() {
    cout << "\nTest SETCHARBITPAIR and GETCHARBITPAIR" << endl;
    unsigned char *array = new unsigned char[1]();
    for (int i = 0; i < 4; i++) {
        cout << " " << GETCHARBITPAIR(array, i);
    }
    cout << endl;
    for (int j = 0; j < 4; j++) {
        for (int k = 1; k < 4; k++) {
            SETCHARBITPAIR(array, j, k);
            for (int i = 0; i < 4; i++) {
                cout << " " << GETCHARBITPAIR(array, i);
            }
            cout << endl;
        }
    }
    for (int j = 0; j < 4; j++) {
        SETCHARBITPAIR(array, j, 0);
        for (int i = 0; i < 4; i++) {
            cout << " " << GETCHARBITPAIR(array, i);
        }
        cout << endl;
    }
}

void testInputData(string inFilename, string outFilename) {
    InputData *input = new InputData(inFilename);
    InputData *output = new InputData(outFilename);
    ULLONG value;
    cout << input->getLength();
    cout << "\n";
    for (ULLONG i = 0; i < input->getLength(); i++) {
        value = input->getValue(i);
        output->setValue(value, i);

        cout << value;
        cout << "\n";
    }
    input->close();
    output->close();
    delete input;
    delete output;
}

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
        cout << "### Verify Perfect Hash Function ###" << endl;
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
            cout << "### Verify Perfect Hash Function Successful ###" << endl;
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
            cout << "### Verify Perfect Hash Function Unsuccessful ###" << endl;
        }
        // Debug end
        return false;
    }
}

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
        << ";size_in_bytes;size_in_bytes_general;size_in_bytes_split_uhf;size_in_bytes_offsets"
        << ";size_in_bytes_good_uhf_pairs;size_in_bytes_random_width;size_in_bytes_random_table"
        << ";size_in_bytes_random_factor;size_in_bytes_g_array"
        <<
        ";compact_size_in_bytes;compact_size_in_bytes_general;compact_size_in_bytes_split_uhf;compact_size_in_bytes_offsets"
        << ";compact_size_in_bytes_good_uhf_pairs;compact_size_in_bytes_random_width;compact_size_in_bytes_random_table"
        << ";compact_size_in_bytes_random_factor;compact_size_in_bytes_g_array"
        << ";creation_time;creation_io"
        << ";setup_time;setup_io"
        << ";split_time;split_io;split_tries;split_success"
        << ";num_of_buckets;max_bucket_size;min_bucket_size;avg_bucket_size"
        << ";goodpairs_time;goodpairs_io;goodpairs_total_tries;goodpairs_success"
        << ";buckets_time;buckets_io;random_tab_tries;random_si_total_tries;buckets_success"
        << ";eval_time;eval_io;avg_eval_time;avg_eval_io;eval_success"
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

    file << ";" << stats.size_in_bytes;
    file << ";" << stats.size_in_bytes_general;
    file << ";" << stats.size_in_bytes_split_uhf;
    file << ";" << stats.size_in_bytes_offsets;
    file << ";" << stats.size_in_bytes_good_uhf_pairs;
    file << ";" << stats.size_in_bytes_random_width;
    file << ";" << stats.size_in_bytes_random_table;
    file << ";" << stats.size_in_bytes_random_factor;
    file << ";" << stats.size_in_bytes_g_array;

    file << ";" << stats.compact_size_in_bytes;
    file << ";" << stats.compact_size_in_bytes_general;
    file << ";" << stats.compact_size_in_bytes_split_uhf;
    file << ";" << stats.compact_size_in_bytes_offsets;
    file << ";" << stats.compact_size_in_bytes_good_uhf_pairs;
    file << ";" << stats.compact_size_in_bytes_random_width;
    file << ";" << stats.compact_size_in_bytes_random_table;
    file << ";" << stats.compact_size_in_bytes_random_factor;
    file << ";" << stats.compact_size_in_bytes_g_array;

    //file << ";" << stats.creation_start;
    //file << ";" << stats.creation_end;
    file << ";" << stats.creation_time;
    file << ";" << stats.creation_io;
    //file << ";" << stats.creation_success;

    //file << ";" << stats.setup_start;
    //file << ";" << stats.setup_end;
    file << ";" << stats.setup_time;
    file << ";" << stats.setup_io;
    //file << ";" << stats.setup_succuess;

    //file << ";" << stats.split_start;
    //file << ";" << stats.split_end;
    file << ";" << stats.split_time;
    file << ";" << stats.split_io;
    file << ";" << stats.split_tries;
    file << ";" << stats.split_success;

    file << ";" << stats.num_of_buckets;
    file << ";" << stats.max_bucket_size;
    file << ";" << stats.min_bucket_size;
    file << ";" << stats.avg_bucket_size;

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

    file << endl;
    file.close();
}

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

    // TODO delete debug code
/*    testInputData();/*
    createInputData(dataFileName);/**/
    createRandomInputData(dataFileName);/**/

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