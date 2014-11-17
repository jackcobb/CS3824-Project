/* 
 * File:   IDnaRepository.h
 * Author: Eric
 *
 * Created on November 3, 2014, 2:07 PM
 */

#ifndef IDNAREPOSITORY_H
#define	IDNAREPOSITORY_H

#include "IDnaSequence.h"
#include "DnaSequenceVector.h"


class IDnaRepository {
public:
    virtual int Size() = 0;
    virtual int Size(int i) = 0;
    virtual int Count() = 0;
    virtual int Count(Nucleotide_t input) = 0;
    virtual int Count(int sequence, Nucleotide_t input) = 0;
    virtual IDnaSequence& Get(int i) =0;
    virtual Nucleotide_t Get(int sequence, int position) =0;
    virtual void Set(int sequence, int position, Nucleotide_t input) = 0;
    virtual void Add(DnaSequenceVector& input) = 0;
private:

};

#endif	/* IDNAREPOSITORY_H */

