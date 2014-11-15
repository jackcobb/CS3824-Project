/* 
 * File:   DnaSequenceVector.h
 * Author: Eric
 *
 * Created on November 15, 2014, 3:03 PM
 */

#ifndef DNASEQUENCEVECTOR_H
#define	DNASEQUENCEVECTOR_H

#include "IDnaSequence.h"


class DnaSequenceVector : IDnaSequence{
public:
    DnaSequenceVector(int maxSize);
    int Size();
    Nucleotide_t Get(int i);
    void Set(int i, Nucleotide_t input);
    private:
        std::vector<Nucleotide_t> members;
};

#endif	/* DNASEQUENCEVECTOR_H */

