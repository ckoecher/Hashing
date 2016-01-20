#include "definitions.h"
#include "PerfectHashFunction.h"

int main() {
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
    return 0;
}