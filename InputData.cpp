/*
 * Hashing - A simple implementation of Split-And-Share-Hashing.
 * Copyright (C) 2016  Philipp Schlag, Chris Köcher
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string.h>
#include "definitions.h"
#include "InputData.h"

using namespace std;

InputData::InputData(string fileName, ios::openmode flags) {
    clock_t time = clock();

    _stream.open(fileName, flags | fstream::in | fstream::out);
    _size = sizeof(ULLONG);

    if(_stream.fail()) {
        cerr << "Could not open file \"" << fileName <<"\". Error: " << errno;
        throw 0;
    }

    //fetch the data length
    _stream.seekg(0, fstream::end);
    _length = (ULLONG) _stream.tellg() / _size;
    _stream.seekg(0);

    _evalTime = clock() - time;
}

InputData::InputData(ios::openmode flags) {
    clock_t time = clock();

    _size = sizeof(ULLONG);
    _length = 0;

    char *tmpName = strdup("/tmp/tmpfileXXXXXX");
    mkstemp(tmpName);
    _stream.open(tmpName, flags | fstream::in | fstream::out);

    if(_stream.fail()) {
        cerr << "Could not open temporary file. Error: " << errno;
        throw 0;
    }

    _evalTime = clock() - time;
}

void InputData::setValue(ULLONG value, ULLONG position) {
    clock_t time = clock();
    char str[_size];

    //convert the data
    for (int i = _size - 1; i >= 0; i--) {
        str[i] = (char) value & 255;
        value >>= 8;
    }

    //write the data
    _stream.seekp((long) position * _size);
    _stream.write(str, _size);

    _evalTime += clock() - time;
}

void InputData::setNextValue(ULLONG value) {
    clock_t time = clock();
    char str[_size];

    //convert the data
    for (int i = _size - 1; i >= 0; i--) {
        str[i] = (char) value & 255;
        value >>= 8;
    }

    //write the data
    _stream.write(str, _size);

    _evalTime += clock() - time;
}

ULLONG InputData::getValue(ULLONG position) {
    clock_t time = clock();
    ULLONG value = 0;
    char str[_size];

    //read the data
    _stream.seekg((long) position * _size);
    _stream.read(str, _size);

    //convert the data
    for (int i = 0; i < _size - 1; i++) {
        value |= (ULLONG) (unsigned char) str[i];
        value <<= 8;
    }
    value |= (ULLONG) (unsigned char) str[_size - 1];

    _evalTime += clock() - time;

    return value;
}

ULLONG InputData::getNextValue() {
    clock_t time = clock();
    ULLONG value = 0;
    char str[_size];

    //read the data
    _stream.read(str, _size);

    //convert the data
    for (int i = 0; i < _size - 1; i++) {
        value |= (ULLONG) (unsigned char) str[i];
        value <<= 8;
    }
    value |= (ULLONG) (unsigned char) str[_size - 1];

    _evalTime += clock() - time;

    return value;
}

ULLONG InputData::getLength() {
    return _length;
}

clock_t InputData::getEvalTime() {
    return _evalTime;
}

void InputData::resetEvalTime() {
    _evalTime = 0;
}

void InputData::close() {
    _stream.close();
}