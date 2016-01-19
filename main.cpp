#include "definitions.h"
#include "PerfectHashFunction.h"

int main() {
    cout << sizeof(unsigned char) << endl;
    unsigned long seed;
    cout << "seed = ";
    cin >> seed;
    mt19937* rng;
    rng = new mt19937(seed);
    uniform_int_distribution<ULLONG>* dist_test;
    dist_test = new uniform_int_distribution<ULLONG>(0,10);
    for(int i=0; i < 100; i++) {
        cout << (*dist_test)(*rng) << endl;
    }
    return 0;
}