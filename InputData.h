//
// Created by chris on 15.02.16.
//
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

#ifndef HASHING_INPUTDATA_H
#define HASHING_INPUTDATA_H

#include "definitions.h"

/**
 * This class is used for handling of the input data.
 */
class InputData {
private:
    /**
     * The internal used file stream.
     */
    std::fstream _stream;

    /**
     * The size of a single item in the data stream.
     */
    int _size;

    /**
     * The number of the items in the data stream.
     */
    ULLONG _length;

    /**
     * The number of ticks used for operations in this class.
     */
    clock_t _evalTime;

public:
    /**
     * Constructor of this class. This method opens the file with the given name.
     *
     * @param[in] fileName The name of the file to open here
     * @param[in] flags    The flags for opening the file
     * @see       ios::openmode
     */
    InputData(std::string fileName, ios::openmode flags = ios::in | ios::out);

    /**
     * Constructor of this class. This method creates a temporary file.
     *
     * @param[in] flags The flags for opening the file
     * @see       ios::openmode
     */
    InputData(ios::openmode flags = ios::in | ios::out);

    /**
     * Inserts an item on the given position in the file.
     *
     * @param[in] value    The value to insert into the file
     * @param[in] position The position in the file
     */
    void setValue(ULLONG value, ULLONG position);

    /**
     * Inserts an item on the position following on the last writing action.
     *
     * @param[in] value The value to insert into the file
     */
    void setNextValue(ULLONG value);

    /**
     * Reads the item on the given position from the file.
     *
     * @param[in]  position The position of the item in file
     * @return              The value of the item
     */
    ULLONG getValue(ULLONG position);

    /**
     * Reads the item on the position following on the last reading action.
     *
     * @param[in]  position The position of the item in file
     * @return              The value of the item
     */
    ULLONG getNextValue();

    /**
     * Returns the number of items that are saved in this file.
     *
     * @return The number of items
     */
    ULLONG getLength();

    /**
     * Returns the number of ticks that are utilized for operations on this file.
     *
     * @return The number of ticks
     */
    clock_t getEvalTime();

    /**
     * Resets the number of ticks.
     */
    void resetEvalTime();

    /**
     * Closes this file. Consider that further operations on this file will raise errors.
     * This method is called automatically on destruction process.
     */
    void close();

    /**
     * Destructor of this class. Just calls the close operation.
     */
    virtual ~InputData() { close(); }
};


#endif //HASHING_INPUTDATA_H
