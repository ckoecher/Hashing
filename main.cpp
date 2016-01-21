#include "definitions.h"
#include "PerfectHashFunction.h"

int main() {
    cout << "Bytes of char: " << sizeof(char) << endl;
    assert(sizeof(char) == 1);
    cout << "Bytes of ULLONG: " << sizeof(ULLONG) << endl;
    assert(sizeof(ULLONG) == 8);
    cout << "Bytes of double: " << sizeof(double) << ", long double: " << sizeof(long double) << endl;
    assert(sizeof(long double) == 16);
    unsigned long seed;
    cout << "seed = ";
    cin >> seed;
    mt19937* rng;
    rng = new mt19937(seed);
    uniform_int_distribution<ULLONG>* dist_test;
    dist_test = new uniform_int_distribution<ULLONG>(0,1);
    for(int i=0; i < 10; i++) {
        cout << (*dist_test)(*rng) << endl;
    }
    cout << "Test of CHARBITPAIR:" << endl;
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
    return 0;
}