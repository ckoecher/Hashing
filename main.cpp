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
    mt19937* rng;
    rng = new mt19937(seed);
    uniform_int_distribution<ULLONG>* dist_test;
    dist_test = new uniform_int_distribution<ULLONG>(0,1);
    for(int i=0; i < 10; i++) {
        cout << (*dist_test)(*rng) << endl;
    }
    delete rng;
    delete dist_test;
}

void testOldBitpairs() {
    cout << "\nTest of CHARBITPAIR:" << endl;
    unsigned char* c = new unsigned char[16];
    cout << "\tchar array:";
    for(unsigned char i = 0; i < 16; i++) {
        c[i] = i;
        cout << " " << i;
    }
    cout << "\n\tCHARBITPAIRs:";
    for(int i = 0; i < 16*4; i++) {
        cout << " " << CHARBITPAIR(c, i);
    }
    cout << "\n\tCHARBITPAIRTABS:";
    for(int i = 0; i < 16*2; i++) {
        cout << " " << CHARBITPAIRTABS(c, 0, i) << " " << CHARBITPAIRTABS(c, 1, i);
    }
    cout << "\nn = ";
    int n;
    clock_t cl;
    unsigned char* ctest;
    unsigned char cc;
    cin >> n;
    ctest = new unsigned char[n];
    for(int i = 0; i < 4*n; i++) {
        cc = CHARBITPAIR3(ctest, i); // caching...
    }
    for(int j = 0; j < 5; j++) {
        cout << "CHARBITPAIR:" << endl;
        cout << "Without optimization: ";
        cl = clock();
        for(int i = 0; i < 4*n; i++) {
            cc = CHARBITPAIR3(ctest, i);
        }
        cout << clock() - cl << " ticks" << endl;
        cout << "With >> 2 instead of / 4: ";
        cl = clock();
        for(int i = 0; i < 4*n; i++) {
            cc = CHARBITPAIR2(ctest, i);
        }
        cout << clock() - cl << " ticks" << endl;
        cout << "With >> 2 instead of / 4 and & 3 instead of % 4: ";
        cl = clock();
        for(int i = 0; i < 4*n; i++) {
            cc = CHARBITPAIR(ctest, i);
        }
        cout << clock() - cl << " ticks" << endl;
        cout << "CHARBITPAIRTABS:" << endl;
        cout << "With >> 2 instead of / 4 and & 3 instead of % 4: ";
        cl = clock();
        for(int i = 0; i < 2*n; i++) {
            cc = CHARBITPAIRTABS(ctest, 0, i);
            cc = CHARBITPAIRTABS(ctest, 1, i);
        }
        cout << clock() - cl << " ticks" << endl;
    }
    unsigned char* charArray;
    cout << "new unsigned char[100]:";
    charArray = new unsigned char[100];
    for(int i = 0; i < 100; i++) {
        cout << " " << (int)charArray[i];
    }
    cout << "\nnew unsigned char[100]():";
    delete[] charArray;
    charArray = new unsigned char[100]();
    for(int i = 0; i < 100; i++) {
        cout << " " << (int)charArray[i];
    }
    int* intArray;
    cout << "\nnew int[100]:";
    intArray = new int[100];
    for(int i = 0; i < 100; i++) {
        cout << " " << intArray[i];
    }
    cout << "\nnew int[100]():";
    delete[] intArray;
    intArray = new int[100]();
    for(int i = 0; i < 100; i++) {
        cout << " " << intArray[i];
    }
    delete[] c;
    delete[] ctest;
    delete[] charArray;
    delete[] intArray;
}

void testBitPairs() {
    cout << "\nNew BitPair implementation: carray[16]";
    int n = 16;
    clock_t cl;
    unsigned char* carray = new unsigned char[n]();
    unsigned char c;
    cout << "carray:";
    for(int i = 0; i < n; i++) {
        carray[i] = i;
        cout << carray[i];
    }
    cout << "\nGETBITPAIR content:";
    for(int i = 0; i < 2*n; i++) {
        for(int j = 0; j < 2; j++) {
            cout << " " << GETBITPAIR(carray, j, i);
        }
        cout << " #";
    }
    cout << "\ngetBitPair content:";
    for(int i = 0; i < 2*n; i++) {
        for(unsigned short j = 0; j < 2; j++) {
            cout << " " << (int)getBitPair(carray, j, i);
        }
        cout << " #";
    }
    cout << "\nZEROBITPAIRS content:";
    for(int i = 0; i < 2*n; i++) {
        for(int j = 0; j < 2; j++) {
            ZEROBITPAIRS(carray, j, i);
        }
    }
    for(int i = 0; i < 2*n; i++) {
        for(int j = 0; j < 2; j++) {
            cout << " " << GETBITPAIR(carray, j, i);
        }
        cout << " #";
    }
    cout << "\nINCBITPAIR content:";
    for(int i = 0; i < 2*n; i++) {
        for(int j = 0; j < 2; j++) {
            if(i % 2 == 1) {
                INCBITPAIR(carray, j, i);
                INCBITPAIR(carray, j, i);
            }
            if(j == 1) {
                INCBITPAIR(carray, j, i);
            }
        }
    }
    for(int i = 0; i < 2*n; i++) {
        for(int j = 0; j < 2; j++) {
            cout << " " << GETBITPAIR(carray, j, i);
        }
        cout << " #";
    }
    cout << "\nn = ";
    cin >> n;
    delete[] carray;
    carray = new unsigned char[n];
    for(int round = 0; round < 5; round++) {
        cout << "\nGETBITPAIR time to read: ";
        cl = clock();
        for(int i = 0; i < 2*n; i++) {
            for(int j = 0; j < 2; j++) {
                c = GETBITPAIR(carray, j, i);
            }
        }
        cout << clock() - cl << " ticks";
        cout << "\ngetBitPair time to read: ";
        cl = clock();
        for(int i = 0; i < 2*n; i++) {
            for(short j = 0; j < 2; j++) {
                c = getBitPair(carray, j, i);
            }
        }
        cout << clock() - cl << " ticks";
    }
    delete[] carray;
}

void testInputData() {
    //InputData *input = new InputData("/home/chris/test.txt");
    //InputData *output = new InputData("/home/chris/test2.txt");
    InputData *input = new InputData("/home/philipp/test.txt");
    InputData *output = new InputData("/home/philipp/test2.txt");
    ULLONG value;
    cout << input->getLength();
    cout << "\n";
    for(ULLONG i = 0; i < input->getLength(); i++) {
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

void testCreateInputData(string dataFileName) {
    cout << "\nCreate Input Data" << endl;
    cout << "n = ";
    ULLONG n;
    cin >> n;
    InputData *data = new InputData(dataFileName, ios::trunc);
    for(ULLONG i = 0; i < n; i++) {
        data->setNextValue(i);
    }
    data->close();
    delete data;
}

void testCreateRandomInputData(string dataFileName) {
    cout << "\nCreate Random Input Data" << endl;
    cout << "n = ";
    ULLONG n;
    cin >> n;
    InputData *data = new InputData(dataFileName, ios::trunc);
    mt19937 *rng = new mt19937(1);
    uniform_int_distribution<ULLONG> *dist = new uniform_int_distribution<ULLONG>;
    for(ULLONG i = 0; i < n; i++) {
        data->setNextValue((*dist)(*rng));
    }
    data->close();
    delete data;
    delete rng;
    delete dist;
}

void testSetGetBit() {
    cout << "\nTest SETBIT and GETBIT" << endl;
    unsigned char *array = new unsigned char[2]();
//#define GETBIT(array, index) ((array[index >> 3]) >> (index & 7)) & 1
//#define SETBIT(array, index, value) array[index >> 3] ^= (-value ^ array[index >> 3]) & (1 << (index & 7))
    for(int i = 0; i < 2*sizeof(unsigned char)*8; i++) {
        cout << GETBIT(array, i) << " ";
    }
    cout << endl;
    for(int j = 0; j < 2*sizeof(unsigned char)*8; j++) {
        SETBIT(array, j, 1);
        for(int i = 0; i < 2*sizeof(unsigned char)*8; i++) {
            cout << GETBIT(array, i) << " ";
        }
        cout << endl;
    }
    for(int j = 0; j < 2*sizeof(unsigned char)*8; j++) {
        SETBIT(array, j, 0);
        for(int i = 0; i < 2*sizeof(unsigned char)*8; i++) {
            cout << GETBIT(array, i) << " ";
        }
        cout << endl;
    }
}

void testSetGetCharBitPairs() {
//#define GETCHARBITPAIR(array, index) ((array[index >> 2]) >> ((index & 3) << 1) & 3)
//#define SETCHARBITPAIR(array, index, value) (array[index >> 2] ^= (-value ^ array[index >> 2]) & (3 << ((index & 3) << 1)))
    cout << "\nTest SETCHARBITPAIR and GETCHARBITPAIR" << endl;
    unsigned char *array = new unsigned char[1]();
    for(int i = 0; i < 4; i++) {
        cout << " " << GETCHARBITPAIR(array, i);
    }
    cout << endl;
    for(int j = 0; j < 4; j++) {
        for(int k = 1; k < 4; k++) {
            SETCHARBITPAIR(array, j, k);
            for(int i = 0; i < 4; i++) {
                cout << " " << GETCHARBITPAIR(array, i);
            }
            cout << endl;
        }
    }
    for(int j = 0; j < 4; j++) {
        SETCHARBITPAIR(array, j, 0);
        for(int i = 0; i < 4; i++) {
            cout << " " << GETCHARBITPAIR(array, i);
        }
        cout << endl;
    }
}

void readConfigs(char* filename, Configuration* configs, ULLONG* numOfConfigs) {
    // reads configuration data from "filename"
    // numOfConfigs configurations -> configs = new Configuration[numOfConfigs]
}

void readData(char* filename, ULLONG* data, ULLONG* data_length) {
    // reads one set of data
}

Configuration readConfig() {
    //INIReader reader("/home/chris/CLionProjects/Hashing/config.ini"); //TODO this path should be changed!!!
    INIReader reader("/home/philipp/ClionProjects/Hashing/config.ini");

    if (reader.ParseError() < 0) {
        cout << "Can't load 'config.ini'\n" << reader.ParseError();
        throw 0; //TODO throw an exception here!
    }

    struct Configuration config;
    config.k = (unsigned short) reader.GetInteger("Hashing", "k", 32);
    config.l = (unsigned short) reader.GetInteger("Hashing", "l", 2);
    config.m_coeff = reader.GetReal("Hashing", "m_coeff", 2);
    config.m_exp = reader.GetReal("Hashing", "m_exp", 2.0/3);
    config.additional_bits_uhf = (unsigned short) reader.GetInteger("Hashing", "additional_bits_uhf", 6);
    config.num_of_tries_split = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_split", 42);
    config.num_of_tries_goodpairs = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_goodpairs", 42);
    config.mi_coeff = reader.GetReal("Hashing", "mi_coeff", 1.25);
    config.tab_rows_coeff = reader.GetReal("Hashing", "tab_rows_coeff", 2);
    config.tab_rows_exp = reader.GetReal("Hashing", "tab_rows_exp", 0.75);
    config.additional_bits_tab = (unsigned short) reader.GetInteger("Hashing", "additional_bits_tab", 6);
    config.num_of_tries_random_tab = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_random_tab", 42); //TODO which default value?
    config.num_of_tries_random_si = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_random_si", 42); //TODO which default value?
    config.seed = (ULLONG) reader.GetInteger("Hashing", "seed", 123456); //TODO which default value?
    config.debug_mode = reader.GetBoolean("Hashing", "debug_mode", true);

    return config;
}

Configuration readConfig(string fileName) {
    INIReader reader(fileName);

    if (reader.ParseError() < 0) {
        //cout << "Can't load 'config.ini'\n" << reader.ParseError();
        cout << "Cannot load configuration file " << fileName << "." << endl << reader.ParseError() << endl;
        throw 0; //TODO throw an exception here!
    }

    struct Configuration config;
    config.k = (unsigned short) reader.GetInteger("Hashing", "k", 32);
    config.l = (unsigned short) reader.GetInteger("Hashing", "l", 2);
    config.m_coeff = reader.GetReal("Hashing", "m_coeff", 2);
    config.m_exp = reader.GetReal("Hashing", "m_exp", 2.0/3);
    config.additional_bits_uhf = (unsigned short) reader.GetInteger("Hashing", "additional_bits_uhf", 6);
    config.num_of_tries_split = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_split", 42);
    config.num_of_tries_goodpairs = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_goodpairs", 42);
    config.mi_coeff = reader.GetReal("Hashing", "mi_coeff", 1.25);
    config.tab_rows_coeff = reader.GetReal("Hashing", "tab_rows_coeff", 2);
    config.tab_rows_exp = reader.GetReal("Hashing", "tab_rows_exp", 0.75);
    config.additional_bits_tab = (unsigned short) reader.GetInteger("Hashing", "additional_bits_tab", 6);
    config.num_of_tries_random_tab = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_random_tab", 42); //TODO which default value?
    config.num_of_tries_random_si = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_random_si", 42); //TODO which default value?
    config.seed = (ULLONG) reader.GetInteger("Hashing", "seed", 123456); //TODO which default value?
    config.debug_mode = reader.GetBoolean("Hashing", "debug_mode", true);

    return config;
}

bool testPerfectHashFunction(Configuration &config, InputData *data, PerfectHashFunction *phf, Statistics &stats) {

    ULLONG datalength = data->getLength();
    ULLONG range = phf->getRange();
    ULLONG key;
    ULLONG hashvalue;
    ULLONG i;
    unsigned char* testarray = new unsigned char[(range>>3)+1]();
    bool perfect = true;

    // Stats
    stats.eval_start = clock();
    // Stats end

    // Debug
    if(config.debug_mode) {
        cout << "### Test Perfect Hash Function ###" << endl;
        cout << "Check " << datalength << " keys:" << endl;
    }
    short percentage = -1;
    // Debug end

    for(i = 0; i < datalength; i++) {
        // Debug
        if(config.debug_mode) {
            if( (i*100) / datalength != percentage) {
                percentage = (i*100) / datalength;
                cout << "\r " << percentage << "% checked." << flush;
            }
        }
        // Debug end
        key = data->getValue(i);
        hashvalue = phf->evaluate(key);
        if(GETBIT(testarray, hashvalue)) {
            perfect = false;
            break;
        }
        SETBIT(testarray, hashvalue, 1);
    }
    // Stats
    stats.eval_end = clock();
    stats.eval_time = stats.eval_end - stats.eval_start;
    // Stats end
    if(perfect) {
        // Debug
        if(config.debug_mode) {
            cout << "\r 100% checked." << endl;
            cout << "### Test Perfect Hash Function Successful ###" << endl;
        }
        // Debug end
        // Stats
        stats.avg_eval_time = (long double)stats.eval_time / (long double)datalength;
        stats.eval_success = true;
        // Stats end
        return true;
    } else {
        // Debug
        if(config.debug_mode) {
            cout << "\nKey no. " << i << " (" << key << ") with  hash value " << hashvalue << " collided." << endl;
            cout << "### Test Perfect Hash Function Unsuccessful ###" << endl;
        }
        // Debug end
        return false;
    }
}

void saveStatistics(string statsFileName, string configFileName, string dataFileName, Statistics stats) {
    ofstream file(statsFileName, ios::out | ios::app);

    file << "\"" << configFileName << "\"";
    file << ";\"" << dataFileName << "\"";

    file << ";" << stats.clocks_per_sec;
    file << ";" << stats.num_of_keys;
    file << ";" << stats.range_of_phf;
    file << ";" << stats.size_in_bytes;

    //file << ";" << stats.creation_start;
    //file << ";" << stats.creation_end;
    file << ";" << stats.creation_time;
    //file << ";" << stats.creation_success;

    //file << ";" << stats.setup_start;
    //file << ";" << stats.setup_end;
    file << ";" << stats.setup_time;
    //file << ";" << stats.setup_succuess;

    //file << ";" << stats.split_start;
    //file << ";" << stats.split_end;
    file << ";" << stats.split_time;
    file << ";" << stats.split_tries;
    file << ";" << stats.split_success;

    file << ";" << stats.num_of_buckets;
    file << ";" << stats.max_bucket_size;
    file << ";" << stats.min_bucket_size;
    file << ";" << stats.avg_bucket_size;

    //file << ";" << stats.goodpairs_start;
    //file << ";" << stats.goodpairs_end;
    file << ";" << stats.goodpairs_time;
    file << ";" << stats.goodpairs_total_tries;
    file << ";" << stats.goodpairs_success;

    //file << ";" << stats.buckets_start;
    //file << ";" << stats.buckets_end;
    file << ";" << stats.buckets_time;
    file << ";" << stats.random_tab_tries;
    file << ";" << stats.random_si_total_tries;
    file << ";" << stats.buckets_success;

    //file << ";" << stats.eval_start;
    //file << ";" << stats.eval_end;
    file << ";" << stats.eval_time;
    file << ";" << stats.avg_eval_time;
    file << ";" << stats.eval_success;

    file << endl;
    file.close();
}

int main(int argc, char* argv[]) {

    // Parse command line arguments
    // argv[0]: path and name of program itself
    // argv[1]: path of configuration file
    // argv[2]: path of data file
    // argv[2]: path of statistics file

//    string configFileName = "/home/philipp/ClionProjects/Hashing/config.ini";
//    //string configFileName = "/home/chris/CLionProjects/Hashing/config.ini";
//    string dataFileName = "/home/philipp/test.txt";
//    //string dataFileName = "/home/chris/test.txt";
//    string statsFileName = "/home/philipp/statistics.csv";
//    //string statsFileName = "/home/chris/statistics.csv";

    string configFileName;
    string dataFileName;
    string statsFileName;

    if(argc < 4) {
        cerr << "Too few command line arguments!" << endl;
        cerr << "./Hashing <configuration file> <data file> <statistics file>" << endl;
        return 1;
    }

    configFileName = argv[1];
    dataFileName = argv[2];
    statsFileName = argv[3];

//    if(argc >= 2) {
//        configFileName = argv[1];
//        if(argc >= 3) {
//            dataFileName = argv[2];
//            if(argc >= 4) {
//                statsFileName = argv[3];
//            }
//        }
//    }

/*    testInputData();/*
    testCreateInputData(dataFileName);/**/
    testCreateRandomInputData(dataFileName);/**/

    // Setup
    Configuration config = readConfig(configFileName);
    InputData *data = new InputData(dataFileName);
    Statistics stats;

    // Creation and Test
    PerfectHashFunction* phf;
    try {
        phf = new PerfectHashFunction(config, data, stats);
        testPerfectHashFunction(config, data, phf, stats);
        delete phf;
    } catch(int e) {
        cerr << "Could not create Perfect Hash Function." << endl;
    }

    // Save statistics
    saveStatistics(statsFileName, configFileName, dataFileName, stats);

    // Cleanup
    delete data;

    return 0;
}