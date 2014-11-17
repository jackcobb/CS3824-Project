/* 
 * File:   DnaSequenceVector.h
 * Author: Eric
 *
 * Created on November 15, 2014, 3:03 PM
 */

#ifndef DNASEQUENCEVECTOR_H
#define	DNASEQUENCEVECTOR_H

#include "IDnaSequence.h"
#include <vector>


class DnaSequenceVector : public IDnaSequence{
public:
    DnaSequenceVector();
    int Size();
    Nucleotide_t Get(int i);
    void Set(int i, Nucleotide_t input);
    void Push(Nucleotide_t input);
    ~DnaSequenceVector();
    private:
       std::vector<Nucleotide_t> members;
};

#endif	/* DNASEQUENCEVECTOR_H */

