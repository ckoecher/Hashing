//
// Created by chris on 15.02.16.
//

#ifndef HASHING_INPUTDATA_H
#define HASHING_INPUTDATA_H

#include "definitions.h"

class InputData {
private:
    std::fstream _stream;
    int _size;
    ULLONG _length;

public:
    InputData(std::string fileName);
    InputData(); //construct a temporary file
    void setValue(ULLONG value, ULLONG position);
    ULLONG getValue(ULLONG position);
    ULLONG getNextValue();
    ULLONG getLength();
    void close();

    virtual ~InputData() { close(); }
};


#endif //HASHING_INPUTDATA_H
