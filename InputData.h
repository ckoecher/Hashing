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
    clock_t _evalTime;

public:
    InputData(std::string fileName, ios::openmode flags = ios::in | ios::out);
    InputData(ios::openmode flags = ios::in | ios::out); //construct a temporary file
    void setValue(ULLONG value, ULLONG position);
    void setNextValue(ULLONG value);
    ULLONG getValue(ULLONG position);
    ULLONG getNextValue();
    ULLONG getLength();
    clock_t getEvalTime();
    void resetEvalTime();
    void close();

    virtual ~InputData() { close(); }
};


#endif //HASHING_INPUTDATA_H
