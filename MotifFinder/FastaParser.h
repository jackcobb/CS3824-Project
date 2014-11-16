/* 
 * File:   FastaParser.h
 * Author: Eric
 *
 * Created on November 3, 2014, 5:48 PM
 */

#ifndef FASTAPARSER_H
#define	FASTAPARSER_H
#include <iostream>     
#include <fstream>

#include "IDnaRepository.h"
class FastaParser {
public:
    FastaParser();
    FastaParser(const FastaParser& orig);
    void Parse(std::ifstream& file, IDnaRepository& repo);
    ~FastaParser();
private:

};

#endif	/* FASTAPARSER_H */

