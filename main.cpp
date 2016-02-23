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

void testCreateInputData() {
    cout << "\nCreate Input Data" << endl;
    cout << "n = ";
    ULLONG n;
    cin >> n;
    //InputData *data = new InputData("/home/chris/test.txt");
    InputData *data = new InputData("/home/philipp/test.txt");
    for(ULLONG i = 0; i < n; i++) {
        data->setNextValue(i);
    }
    data->close();
    delete data;
}

void testCreateRandomInputData() {
    cout << "\nCreate Random Input Data" << endl;
    cout << "n = ";
    ULLONG n;
    cin >> n;
    //InputData *data = new InputData("/home/chris/test.txt");
    InputData *data = new InputData("/home/philipp/test.txt");
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
    config.mi_coeff = reader.GetReal("Hashing", "mi_coeff", 1.25);
    config.tab_rows_coeff = reader.GetReal("Hashing", "tab_rows_coeff", 2);
    config.tab_rows_exp = reader.GetReal("Hashing", "tab_rows_exp", 0.75);
    config.additional_bits_tab = (unsigned short) reader.GetInteger("Hashing", "additional_bits_tab", 6);
    config.num_of_tries_random_tab = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_random_tab", 42); //TODO which default value?
    config.num_of_tries_random_si = (unsigned short) reader.GetInteger("Hashing", "num_of_tries_random_si", 42); //TODO which default value?
    config.seed = (ULLONG) reader.GetInteger("Hashing", "seed", 123456); //TODO which default value?

    return config;
}

int main() {
/*    testBitPairs();
    testTypeSizes();
    testMersenneTwister();
    testOldBitpairs();*/

    Configuration config = readConfig();
/*    ULLONG data_length = 10;
    ULLONG* data = new ULLONG[data_length];
    for(ULLONG i = 0; i < data_length; i++) {
        data[i] = i;
    }*/


    // TODO n = 16 => "terminate called after throwing an instance of 'int'" after createGoodPairs
    // TODO exists bucket with >= 2 elements that is acyclic (by function)?
    // TODO n = 5000 => bucket 0 no good pair of hash functions
/*    testInputData();/**/
    testCreateInputData();/*
    testCreateRandomInputData();/**/
    //InputData *data = new InputData("/home/chris/test.txt");
    InputData *data = new InputData("/home/philipp/test.txt");

    PerfectHashFunction* phf = new PerfectHashFunction(config, data);

    data->close();
    delete data;
    delete phf;

    return 0;
}