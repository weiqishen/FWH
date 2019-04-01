/*! \class param_reader
 *  \brief Simple, robust method for reading input files
 *  \author Jacob Crabill
 *  \date 4/30/2015
 */
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include "ndarray.h"

class param_reader
{
public:

    /*! Default constructor */
    param_reader();

    /**
     * @brief Construct a new param reader object
     * 
     * @param input file to read from 
     */
    param_reader(string fileName);

    /*! Default destructor */
    ~param_reader();

    /*! Set the file to be read from */
    void setFile(string fileName);

    /*! Open the file to prepare for reading simulation parameters */
    void openFile(void);

    /*! Close the file & clean up */
    void closeFile(void);

    /* === Functions to read paramters from input file === */

    /*! Read a single value from the input file; if not found, apply a default value */
    template <typename T>
    void getScalarValue(string optName, T &opt, T defaultVal);

    /*! Read a single value from the input file; if not found, throw an error and exit */
    template <typename T>
    void getScalarValue(string optName, T &opt);

    /*! Read a vector of values from the input file; if not found, apply the default value to all elements */
    template <typename T>
    void getVectorValue(string optName, vector<T> &opt, T defaultVal);

    /*! Read a vector of values from the input file; if not found, throw an error and exit */
    template <typename T>
    void getVectorValue(string optName, vector<T> &opt);

    /*! Read a ndarray of values from the input file; if not found, throw an error and exit */
    template <typename T>
    void getVectorValue(string optName, ndarray<T> &opt);
private:
    ifstream optFile;
    string fileName;

};


using namespace std;

param_reader::param_reader()
{
}

param_reader::param_reader(string fileName)
{
    this->fileName = fileName;
}

param_reader::~param_reader()
{
    if (optFile.is_open()) optFile.close();
}

void param_reader::setFile(string fileName)
{
    this->fileName = fileName;
}

void param_reader::openFile(void)
{
    optFile.open(fileName.c_str(), ifstream::in);
}

void param_reader::closeFile()
{
    optFile.close();
}

template<typename T>
void param_reader::getScalarValue(string optName, T &opt, T defaultVal)
{
    string str, optKey;

    if (!optFile.is_open() || !getline(optFile,str))
    {
        optFile.open(fileName.c_str());
        if (!optFile.is_open())
            Fatal_Error("Cannont open input file for reading.");
    }

    // Rewind to the start of the file
    optFile.clear();
    optFile.seekg(0,optFile.beg);

    // Search for the given option string
    while (getline(optFile,str))
    {
        // Remove any leading whitespace & see if first word is the input option
        stringstream ss;
        ss.str(str);
        ss >> optKey;
        if (optKey.compare(optName)==0)
        {
            if (!(ss >> opt))
            {
                // This could happen if, for example, trying to assign a string to a double
                cout << "WARNING: Unable to assign value to option " << optName << endl;
                cout << "Using default value of " << defaultVal << " instead." << endl;
                opt = defaultVal;
            }

            return;
        }
    }

    opt = defaultVal;
}

template<typename T>
void param_reader::getScalarValue(string optName, T &opt)
{
    string str, optKey;

    if (!optFile.is_open())
    {
        optFile.open(fileName.c_str());
        if (!optFile.is_open())
            Fatal_Error("Cannont open input file for reading.");
    }

    // Rewind to the start of the file
    optFile.clear();
    optFile.seekg(0,optFile.beg);

    // Search for the given option string
    while (getline(optFile,str))
    {
        // Remove any leading whitespace & see if first word is the input option
        stringstream ss;
        ss.str(str);
        ss >> optKey;

        if (optKey.compare(optName)==0)
        {
            if (!(ss >> opt))
            {
                // This could happen if, for example, trying to assign a string to a double
                cerr << "WARNING: Unable to assign value to option " << optName << endl;
                string errMsg = "Required option not set: " + optName;
                Fatal_Error(errMsg.c_str())
            }

            return;
        }
    }

    // Option was not found; throw error & exit
    string errMsg = "Required option not found: " + optName;
    Fatal_Error(errMsg.c_str())
}

template<typename T>
void param_reader::getVectorValue(string optName, vector<T> &opt)
{
    string str, optKey;

    if (!optFile.is_open())
    {
        optFile.open(fileName.c_str());
        if (!optFile.is_open())
            Fatal_Error("Cannont open input file for reading.");
    }

    // Rewind to the start of the file
    optFile.clear();
    optFile.seekg(0,optFile.beg);

    // Search for the given option string
    while (getline(optFile,str))
    {
        // Remove any leading whitespace & see if first word is the input option
        stringstream ss;
        ss.str(str);
        ss >> optKey;
        if (optKey.compare(optName)==0)
        {
            int nVals;
            if (!(ss >> nVals))
            {
                // This could happen if, for example, trying to assign a string to a double
                cerr << "WARNING: Unable to read number of entries for vector option " << optName << endl;
                string errMsg = "Required option not set: " + optName;
                Fatal_Error(errMsg.c_str());
            }

            opt.resize(nVals);
            for (int i=0; i<nVals; i++)
            {
                if (!ss >> opt[i])
                {
                    cerr << "WARNING: Unable to assign all values to vector option " << optName << endl;
                    string errMsg = "Required option not set: " + optName;
                    Fatal_Error(errMsg.c_str())
                }
            }

            return;
        }
    }

    // Option was not found; throw error & exit
    string errMsg = "Required option not found: " + optName;
    Fatal_Error(errMsg.c_str())
}

template<typename T>
void param_reader::getVectorValue(string optName, ndarray<T> &opt)
{
    string str, optKey;

    if (!optFile.is_open())
    {
        optFile.open(fileName.c_str());
        if (!optFile.is_open())
            Fatal_Error("Cannont open input file for reading.");
    }

    // Rewind to the start of the file
    optFile.clear();
    optFile.seekg(0,optFile.beg);

    // Search for the given option string
    while (getline(optFile,str))
    {
        // Remove any leading whitespace & see if first word is the input option
        stringstream ss;
        ss.str(str);
        ss >> optKey;
        if (optKey.compare(optName)==0)
        {
            int nVals;
            if (!(ss >> nVals))
            {
                // This could happen if, for example, trying to assign a string to a double
                cerr << "WARNING: Unable to read number of entries for vector option " << optName << endl;
                string errMsg = "Required option not set: " + optName;
                Fatal_Error(errMsg.c_str());
            }

            opt.setup(nVals);
            for (int i=0; i<nVals; i++)
            {
                if (!(ss >> opt(i)))
                {
                    cerr << "WARNING: Unable to assign all values to vector option " << optName << endl;
                    string errMsg = "Required option not set: " + optName;
                    Fatal_Error(errMsg.c_str());
                }
            }

            return;
        }
    }

    // Option was not found; throw error & exit
    string errMsg = "Required option not found: " + optName;
    Fatal_Error(errMsg.c_str())
}
