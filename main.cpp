#include "definitions.h"
#include "PerfectHashFunction.h"

void testTypeSizes() {
    cout << "\nBytes of char: " << sizeof(char) << endl;
    assert(sizeof(char) == 1);
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
    char* c = new char[16];
    cout << "\tchar array:";
    for(char i = 0; i < 16; i++) {
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
    char* ctest;
    char cc;
    cin >> n;
    ctest = new char[n];
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
    char* charArray;
    cout << "new char[100]:";
    charArray = new char[100];
    for(int i = 0; i < 100; i++) {
        cout << " " << (int)charArray[i];
    }
    cout << "\nnew char[100]():";
    delete[] charArray;
    charArray = new char[100]();
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
    char* carray = new char[n]();
    char c;
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
        for(short j = 0; j < 2; j++) {
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
    carray = new char[n];
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

void readConfigs(char* filename, Configuration* configs, ULLONG* numOfConfigs) {
    // reads configuration data from "filename"
    // numOfConfigs configurations -> configs = new Configuration[numOfConfigs]
}

void readData(char* filename, ULLONG* data, ULLONG* data_length) {
    // reads one set of data
}

Configuration readConfig() {
    INIReader reader("/home/chris/CLionProjects/Hashing/config.ini"); //TODO this path should be changed!!!

    if (reader.ParseError() < 0) {
        cout << "Can't load 'config.ini'\n" << reader.ParseError();
        throw 0; //TODO throw an exception here!
    }

    struct Configuration config;
    config.k = (short) reader.GetInteger("Hashing", "k", 8);
    config.l = (short) reader.GetInteger("Hashing", "l", 8);
    config.m_coeff = reader.GetReal("Hashing", "m_coeff", 2);
    config.m_exp = reader.GetReal("Hashing", "m_exp", 2/3);
    config.additional_bits_uhf = (short) reader.GetInteger("Hashing", "additional_bits_uhf", 6);
    config.mi_coeff = reader.GetReal("Hashing", "mi_coeff", 1.25);
    config.tab_rows_coeff = reader.GetReal("Hashing", "tab_rows_coeff", 2);
    config.tab_rows_exp = reader.GetReal("Hashing", "tab_rows_exp", 0.75);
    config.additional_bits_tab = (short) reader.GetInteger("Hashing", "additional_bits_tab", 6);
    config.num_of_tries_random_tab = (short) reader.GetInteger("Hashing", "num_of_tries_random_tab", 42); //TODO which default value?
    config.num_of_tries_random_si = (short) reader.GetInteger("Hashing", "num_of_tries_random_si", 42); //TODO which default value?
    config.seed = (unsigned long) reader.GetInteger("Hashing", "seed", 123456); //TODO which default value?

    return config;
}

int main() {
    readConfig();
/*    testBitPairs();
    testTypeSizes();
    testMersenneTwister();
    testOldBitpairs();*/

    return 0;
}