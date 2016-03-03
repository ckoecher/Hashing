//
// Created by chris on 15.02.16.
//

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
     * The amount of the items in the data stream.
     */
    ULLONG _length;

    /**
     * The count of ticks used for operations in this class.
     */
    clock_t _evalTime;

public:
    /**
     * Constructor of this class. This method opens the file with the given name.
     *
     * @param fileName The name of the file to open here
     * @param flags    The flags for opening the file
     * @see   ios::openmode
     */
    InputData(std::string fileName, ios::openmode flags = ios::in | ios::out);

    /**
     * Constructor of this class. This method creates a temporary file.
     *
     * @param flags The flags for opening the file
     * @see   ios::openmode
     */
    InputData(ios::openmode flags = ios::in | ios::out); //construct a temporary file

    /**
     * Inserts an item on the given position in the file.
     *
     * @param value    The value to insert into the file
     * @param position The position in the file
     */
    void setValue(ULLONG value, ULLONG position);

    /**
     * Inserts an item on the position following on the last writing action.
     *
     * @param value The value to insert into the file
     */
    void setNextValue(ULLONG value);

    /**
     * Reads the item on the given position from the file.
     *
     * @param  position The position of the item in file
     * @return          The value of the item
     */
    ULLONG getValue(ULLONG position);

    /**
     * Reads the item on the position following on the last reading action.
     *
     * @param  position The position of the item in file
     * @return          The value of the item
     */
    ULLONG getNextValue();

    /**
     * Returns the amount of items that are saved in this file.
     *
     * @return The amount of items
     */
    ULLONG getLength();

    /**
     * Returns the amount of ticks that are utilized for operations on this file.
     *
     * @return The amount of ticks
     */
    clock_t getEvalTime();

    /**
     * Resets the amount of ticks.
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
