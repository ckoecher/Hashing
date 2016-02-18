//
// Created by chris on 15.02.16.
//

#include <string.h>
#include "definitions.h"
#include "InputData.h"

using namespace std;

InputData::InputData(string fileName) {
    _stream.open(fileName, fstream::in | fstream::out);
    _size = sizeof(ULLONG);

    //fetch the data length
    _stream.seekg(0, fstream::end);
    _length = (ULLONG) _stream.tellg() / _size;

    _stream.seekg(0);
}

InputData::InputData() {
    _size = sizeof(ULLONG);
    _length = 0;

    char *tmpName = strdup("/tmp/tmpfileXXXXXX");
    mkstemp(tmpName);
    _stream.open(tmpName, fstream::in | fstream::out);
}

void InputData::setValue(ULLONG value, ULLONG position) {
    char str[_size];

    //convert the data
    for(int i = _size - 1; i >= 0; i--) {
        str[i] = (char) value & 255;
        value >>= 8;
    }

    //write the data
    _stream.seekp((long) position * _size);
    _stream.write(str, _size);
}

void InputData::setNextValue(ULLONG value) {
    char str[_size];

    //convert the data
    for(int i = _size - 1; i >= 0; i--) {
        str[i] = (char) value & 255;
        value >>= 8;
    }

    //write the data
    _stream.write(str, _size);
}

ULLONG InputData::getValue(ULLONG position) {
    ULLONG value = 0;
    char str[_size];

    //read the data
    _stream.seekg((long) position * _size);
    _stream.read(str, _size);

    //convert the data
    for(int i = 0; i < _size - 1; i++) {
        value |= (ULLONG) str[i];
        value <<= 8;
    }
    value |= (ULLONG) str[_size - 1];

    return value;
}

ULLONG InputData::getNextValue() {
    ULLONG value = 0;
    char str[_size];

    //read the data
    _stream.read(str, _size);

    //convert the data
    for(int i = 0; i < _size - 1; i++) {
        value |= (ULLONG) str[i];
        value <<= 8;
    }
    value |= (ULLONG) str[_size - 1];

    return value;
}

ULLONG InputData::getLength() {
    return _length;
}

void InputData::close() {
    _stream.close();
}